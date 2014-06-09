========
Usage
========

To use Absolute_Integrator in a project::

	import Absolute_Integrator
        
Absolute Integrator provides several use cases.  Each of these
generally provide a simple function that wraps multiple methods of
achieving a goal.  For example, the peak_finding module::

        import Absolute_Integrator.peak_finding as pf
        pf.list_methods()
        pf.list_options(method="ranger")
        pf.peak_find(data, method="ranger")

