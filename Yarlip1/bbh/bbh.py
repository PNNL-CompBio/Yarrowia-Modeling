import click
from sultan.api import Sultan
import pandas as pd
import os

@click.command()
@click.option('--fasta1', help='First protein fasta')
@click.option('--fasta2', help='Second protein fasta')
@click.option('--outdir', help='Output directory')
def bbh(fasta1, fasta2, outdir):
    "Run Bidirectional best hits on two organisms"
    dir1, dir2 = os.path.dirname(fasta1), os.path.dirname(fasta2)
    base1, base2 = os.path.basename(fasta1), os.path.basename(fasta2)
    org1, suffix1 = os.path.splitext(base1)
    org2, suffix2 = os.path.splitext(base2)
    exonerate_args = """\
        --query {} \
        --target {} \
        --bestn 1 \
        --ryo "%qi\t%ti\t%ps\n" \
        --showvulgar no \
        --verbose 0 \
        --showalignment no\
        """
    outfile = '{}.to.{}.tab'
    header1 = [org1, org2, 'Similarity_{}.to.{}'.format(org1, org2)]
    header2 = [org2, org1, 'Similarity_{}.to.{}'.format(org2, org1)]
    with Sultan.load() as s:
        s.exonerate( exonerate_args.format( fasta1,
                                            fasta2 ) ).redirect(
            os.path.join( outdir,
                          outfile.format( org1,
                                          org2 ) ), 
            append=False,
            stdout=True,
            stderr=False ).run()
        s.exonerate( exonerate_args.format( fasta2,
                                            fasta1 ) ).redirect(
            os.path.join( outdir,
                          outfile.format( org2,
                                          org1 ) ), 
            append=False,
            stdout=True,
            stderr=False).run()
    org1_to_org2 = pd.read_table( os.path.join( outdir,
                                                outfile.format( org1,
                                                                org2 ) ),
                                  names=header1 )
    org2_to_org1 = pd.read_table( os.path.join( outdir,
                                                outfile.format( org2,
                                                                org1 ) ),
                                  names=header2 )
    bbh = org1_to_org2.merge( org2_to_org1,
                              on   = [org1, org2],
                              how  = 'inner',
                              sort = True ).sort_values( by=header1[2],
                                                         ascending=False )
    bbh.to_csv( os.path.join( outdir,
                              '{}_and_{}_BBH.tab'.format( org1,
                                                          org2 ) ),
                index=False,
                sep='\t' )
                              
            
if __name__ == '__main__':
    bbh()
    
