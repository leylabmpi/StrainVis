import panel as pn
from SynTrackerVis_app.syntracker_vis_app import SynTrackerVisApp

MAX_SIZE_MB = 500


def destroyed(session_context):
    print("\n\n\nThe session is closed from syntracker_vis.py...")


def get_user_page():
    app = SynTrackerVisApp()
    template = app.template
    pn.state.on_session_destroyed(destroyed)
    return template


def serve_user_page():
    pn.serve(
        get_user_page(),
        port=5006,
        show=True,
        # Increase the maximum websocket message size allowed by Bokeh
        websocket_max_message_size=MAX_SIZE_MB*1024*1024,
        # Increase the maximum buffer size allowed by Tornado
        http_server_kwargs={'max_buffer_size': MAX_SIZE_MB*1024*1024, 'max_body_size': MAX_SIZE_MB*1024*1024}
    )


# Running script with 'python main.py' - use pn.serve
if __name__ == '__main__':
    serve_user_page()

# Running script with 'panel serve'
else:
    get_user_page()


