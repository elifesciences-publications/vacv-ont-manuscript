from . import print_summary
from . import condensed
from . import hist
from . import sb
from . import __version__

import sys

commands = [
           ('summary', 'text summary of ONT data', print_summary.main),
           ('condensed', 'linked boxes representing allele frequencies within copies', condensed.main),
           ('hist', 'stacked bars representing copy number distributions', hist.main),
           ('sb', 'stacked bars representing frequencies of allele combinations', sb.main),
           ]

def print_commands():
    print "vacv-plot version: %s" % (__version__)
    print("""
commands:
%s""" % "\n".join(("%-14s: %s" % (c[0], c[1])) for c in commands))

def main():
    if len(sys.argv) < 2 or not sys.argv[1] in (c[0] for c in commands):
        print_commands() 
        sys.exit(1)
    
    cmd = next(c for c in commands if c[0] == sys.argv[1])
    cmd[2](sys.argv[2:])

if __name__ == "__main__":
    main() 
