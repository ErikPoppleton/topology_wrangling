"""The one exception type the library raises."""


class TopologyError(Exception):
    """Something is wrong with the input that the user needs to fix.

    Library code never calls sys.exit(); it raises this and lets the command
    line front-end in topology_wrangling.cli.common turn it into a message and
    an exit status.  That is what makes every module here importable from a
    script that wants to handle the problem itself.
    """
