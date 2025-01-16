import panel as pn
import time
import random
from syntracker_vis_app import SynTrackerVisApp


app_instance = SynTrackerVisApp()
app_instance.get_page()
#pn.serve(app_instance.app, port=5000, show=False)


