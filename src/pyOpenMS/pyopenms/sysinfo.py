import sys

if sys.platform.startswith("linux"):
    import ctypes as c
    import ctypes.util

    libc = c.CDLL(ctypes.util.find_library("c"))

    class SysInfo(c.Structure):

         # for libc5: needs padding
        padding = 20 - 2 * c.sizeof(c.c_long) - c.sizeof(c.c_int)

        _fields_ = [("uptime", c.c_ulong),
                    ("loads", 3 * c.c_ulong),
                    ("totalram", c.c_ulong),
                    ("freeram", c.c_ulong),
                    ("sharedram", c.c_ulong),
                    ("bufferram", c.c_ulong),
                    ("totalswap", c.c_ulong),
                    ("freeswap", c.c_ulong),
                    ("procs", c.c_ushort),
                    ("totalhigh", c.c_ulong),
                    ("freehigh", c.c_ulong),
                    ("mem_unit", c.c_ulong),
                    ("_padding", padding * c.c_char)
        ]

    def free_mem():
        sys_info = SysInfo()
        libc.sysinfo(c.byref(sys_info))
        return sys_info.freeram

elif sys.platform == "win32":
    try:
        import win32api
    except:
        free_mem = lambda: 0 # memory will never change !
    else:
        def free_mem():
            return win32api.GlobalMemoryStatus()['AvailPhys']

else:
    sys.stderr.write("Determination of memory status is not supported on this \n"
                     " platform, measuring for memoryleaks will never fail\n")

    free_mem = lambda: 0 # memory will never change !
