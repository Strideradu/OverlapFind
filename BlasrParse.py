class BlasrParse(object):
    def __init__(self, blasr_line):
        line = blasr_line.rstrip()
        line_sp = line.split()

        if len(line_sp) == 13:
            # m4 case

            self.id = line_sp[0]

            self.target_id = line_sp[1]
            self.score = int(line_sp[2])
            self.identity = float(line_sp[3])
            self.query_strand = int(line_sp[4])
            self.query_start = int(line_sp[5])
            self.query_end = int(line_sp[6])
            self.query_len = int(line_sp[7])
            self.target_strand = int(line_sp[8])

            self.target_len = int(line_sp[11])
            if self.target_strand == 0:
                self.target_start = int(line_sp[9])
                self.target_end = int(line_sp[10])
            else:
                self.target_start = self.target_len - int(line_sp[10])
                self.target_end = self.target_len - int(line_sp[9])

            self.QV = int(line_sp[12])