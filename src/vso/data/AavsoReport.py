
class AavsoReport:
    """Report generator

        See https://www.aavso.org/aavso-extended-file-format for format specification
    """
    def __init__(self, file, target, chart, obscode) -> None:
        self.file_ = file
        self.target_ = target.upper()
        self.chart_ = chart
        self.obscode_ = obscode

    def header(self):
       print("#TYPE=Extended", file=self.file_)
       print(f"#OBSCODE={self.obscode_}", file=self.file_)
       print("#SOFTWARE=Astropy-based framework: https://github.com/dmitrymu/vsopy", file=self.file_)
       print("#DELIM=,", file=self.file_)
       print("#DATE=JD", file=self.file_)
       print("#OBSTYPE=CCD", file=self.file_)
       print("#NAME,DATE,"
             "MAG,MERR,"
             "FILT,TRANS,MTYPE,CNAME,CMAG,"
             "KNAME,KMAG,"
             "AMASS,GROUP,CHART,NOTES", file=self.file_)

    def body(self, data, band):
      for obs in data:
       print(f"{self.target_},{obs['time']:14f},"
             f"{obs[band]['mag'].value:6.3f},{obs[band]['err'].value:.2g},"
             f"{band[0]},NO,STD,{obs['comp']},{
                 obs[f'comp {band}']['mag']:6.3f},"
             f"{obs['check']},{obs[f'check {band}']['mag'].value:6.3f},"
             f"{obs['airmass']:.5g},{obs['batch']},{
                 self.chart_},Transform_method=simple",
             file=self.file_)
