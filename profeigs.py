import eigs
import line_profiler

profile = line_profiler.LineProfiler(eigs.main)

profile.runcall(eigs.main)
profile.print_stats()

