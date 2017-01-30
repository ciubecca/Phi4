import builtins
from memory_profiler import memory_usage
import warnings

try:
        profile = builtins.profile
except AttributeError:
        # No line profiler, provide a pass-through version
        def profile(func):
            msg = "Leaving profile decorators for {} might"\
                    " affect performance".format(func)
            warnings.warn(msg, RuntimeWarning)
            return func

