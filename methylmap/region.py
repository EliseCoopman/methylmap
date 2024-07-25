import sys


class Region(object):
    def __init__(self, region=None, tup=None):
        if tup:
            self.chromosome, self.begin, self.end = tup
        else:
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
        self.start = self.begin  # start as alias for begin
        self.size = self.end - self.begin
        if self.size < 0:
            sys.exit(
                "\n\nERROR: Window (-w/--window) inproperly formatted, "
                "begin of the interval has to be smaller than end\n\n"
            )
        self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        self.fmt = f"{self.chromosome}:{self.begin}-{self.end}"

    def __mul__(self, other):
        new_half_size = int((self.size * other) / 2)
        middlepoint = (self.begin + self.end) / 2
        return Region(
            tup=(
                self.chromosome,
                int(middlepoint - new_half_size),
                int(middlepoint + new_half_size),
            )
        )

    def __truediv__(self, other):
        new_half_size = int((self.size / other) / 2)
        middlepoint = (self.begin + self.end) / 2
        return Region(
            tup=(
                self.chromosome,
                int(middlepoint - new_half_size),
                int(middlepoint + new_half_size),
            )
        )

