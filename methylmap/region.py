import sys

class Region(object):
    def __init__(self, region, expand=False):
        try:
            self.chromosome, interval = region.replace(",", "").split(":")
            try:
                # see if just integer chromosomes are used
                self.chromosome = int(self.chromosome)
            except ValueError:
                pass
            self.begin, self.end = [int(i) for i in interval.split("-")]
        except ValueError:
            sys.exit(
                "\n\nERROR: Window (-w/--window) inproperly formatted, "
                "an example of accepted formats is:\n'chr5:150200605-150423790'\n\n"
            )
        if expand:
            self.begin = self.begin - int(expand)
            self.end = self.end + int(expand)
        self.start = self.begin  # start as alias for begin
        self.size = self.end - self.begin
        if self.size < 0:
            sys.exit(
                "\n\nERROR: Window (-w/--window) inproperly formatted, "
                "begin of the interval has to be smaller than end\n\n"
            )
        self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        self.fmt = f"{self.chromosome}:{self.begin}-{self.end}"