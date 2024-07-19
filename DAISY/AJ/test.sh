time awk '
NR==FNR {patterns[$0]; next}
{
  header = $0
  getline seq
  for (pattern in patterns) {
    if (match(seq, pattern) > 0) {
      print header
      print seq
      break
    }
  }
}
' barcode_test.txt SB28_DAISY_DAY0_sub2_A3_S122_L007_unpaired_R1_homdomain.fasta > SB28_DAISY_DAY0_sub2_A3_S122_L007_unpaired_R1_homdomain_10x_barcode_test.fasta