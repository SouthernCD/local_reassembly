class ChrLoci(object):
    def __init__(self, chr_id=None, strand=None, start=None, end=None, sp_id=None):
        self.chr_id = chr_id
        self.sp_id = sp_id

        if strand is None:
            self.strand = strand
        elif strand == "+" or str(strand) == '1':
            self.strand = "+"
        elif strand == "-" or str(strand) == '-1':
            self.strand = "-"
        else:
            self.strand = None

        if end is not None and start is not None:
            self.start = min(int(start), int(end))
            self.end = max(int(start), int(end))
            self._range = (self.start, self.end)
            self.length = abs(self.end - self.start) + 1

    # def __eq__(self, other):
    #     """Implement equality by comparing all the location attributes."""
    #     if not isinstance(other, ChrLoci):
    #         return False
    #     return self.start == other.start and \
    #            self.end == other.end and \
    #            self.strand == other.strand and \
    #            self.chr_id == other.chr_id

    def get_fancy_name(self):
        if self.strand is None:
            self.fancy_name = "%s:%d-%d" % (self.chr_id, self.start, self.end)
        else:
            self.fancy_name = "%s:%d-%d:%s" % (self.chr_id,
                                               self.start, self.end, self.strand)
        return self.fancy_name

    def __str__(self):
        try:
            self.get_fancy_name()
            if self.sp_id:
                return "%s %s" % (self.sp_id, self.fancy_name)
            else:
                return self.fancy_name
        except:
            return "No detail range"

    __repr__ = __str__

    def len(self):
        return abs(self.end - self.start) + 1

    def __eq__(self, other):
        return self.chr_id == other.chr_id and self.strand == other.strand and self.start == other.start and self.end == other.end and self.sp_id == other.sp_id

    def __hash__(self):
        return hash(id(self))

