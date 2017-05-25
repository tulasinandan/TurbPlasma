#!/usr/bin/env python
from subs import create_object
rc=create_object()

print
print rc.dirname
print
print
rc.print_params()
rc.loadenergies()
print 'E_b0=',rc.eb0
