import argparse
import importlib
import time
import panel as pn

RESTART_DELAY = 2  # seconds
MAX_RETRIES = 3
MAX_FILE_SIZE_BYTES = 524288000  # 500 Mb


def start_server(app_filename, port, show):
    """
    Imports the app script dynamically and launches the Panel server
    with custom Tornado HTTP socket buffer sizes.
    """
    # Remove .py extension to convert file name to module path
    module_name = app_filename.replace('.py', '')

    # Import the target module dynamically
    app_module = importlib.import_module(module_name)

    # Retrieve the create_app factory function from the module
    if not hasattr(app_module, 'create_app'):
        raise AttributeError(f"Module '{module_name}' does not define a 'create_app()' function.")

    app_factory = getattr(app_module, 'create_app')

    print(f"\nStarting Panel server for '{app_filename}' on port {port}...")

    pn.serve(
        {'strain_vis_debug': app_factory},
        port=int(port),
        show=show,
        websocket_max_message_size=MAX_FILE_SIZE_BYTES,
        unused_session_lifetime=360000,
        http_server_kwargs={
            'max_buffer_size': MAX_FILE_SIZE_BYTES  # Fixes Firefox read buffer limit
        }
    )


def run_with_retry(app_filename, port, show):
    retry_count = 0
    while True:
        try:
            start_server(app_filename, port, show)
            # If server closes cleanly
            print("\nPanel server stopped. Restarting...")
            time.sleep(RESTART_DELAY)
            retry_count = 0
        except Exception as e:
            retry_count += 1
            print(f"\nPanel crash detected: {e}")
            print(f"Retry {retry_count}/{MAX_RETRIES}")

            if retry_count >= MAX_RETRIES:
                print("\nToo many crashes, giving up...")
                break

            time.sleep(RESTART_DELAY)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=str, default="5005")
    parser.add_argument("--show", action='store_true', default=False)
    args = parser.parse_args()

    app = "strain_vis_debug.py"

    run_with_retry(app, args.port, args.show)
