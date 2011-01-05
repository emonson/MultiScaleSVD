from bbfreeze import Freezer
f = Freezer("mssvd_gui_color_0104", includes=("sip",))
f.addScript("main.py")
f()
