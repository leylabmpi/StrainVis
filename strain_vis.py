import panel as pn
import os
import psutil
from StrainVis_app.strain_vis_app import StrainVisApp

MAX_MB = 2 * 1024  # 2 Gb

# --- Module-level flag to ensure single registration ---
_restart_monitor_registered = False


def get_memory_mb():
    proc = psutil.Process(os.getpid())
    try:
        mem = proc.memory_full_info().uss / 1024 / 1024  # MB
    except AttributeError:
        mem = proc.memory_info().rss / 1024 / 1024  # fallback MB
    return mem


def maybe_restart():
    mem_mb = get_memory_mb()

    # Count active sessions safely
    active_sessions = 0
    for s in pn.state.session_info.values():
        if isinstance(s, dict):
            if s.get("started") and not s.get("ended"):
                active_sessions += 1
        else:
            # fallback: assume it's an active session
            active_sessions += 1

    #print(f"[Monitor] Memory: {mem_mb: .1f} MB | Active sessions: {active_sessions}")

    if mem_mb > MAX_MB and active_sessions <= 1:
        print("\n\nRestarting panel server due to memory limit...")
        os._exit(0)  # supervisor will restart


def register_monitor():
    """Register the memory monitor once per server."""
    global _restart_monitor_registered
    if not _restart_monitor_registered:
        #print("\nRegistering maybe_restart periodic callback on server")
        pn.state.add_periodic_callback(maybe_restart, period=60000)  # every 60 sec
        _restart_monitor_registered = True


def create_app():
    print("\n\nCreating new session, ID:", pn.state.curdoc.session_context.id)

    app = StrainVisApp()

    # Build the UI
    app._build_template()

    # bind cleanup explicitly
    pn.state.on_session_destroyed(app._cleanup)

    # Register monitor once per server process
    register_monitor()

    return app.template


# Create the app for one session
create_app().servable("strain_vis")

