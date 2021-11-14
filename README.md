# Homework 2

EXERCISE 1 DNA
         
          def countNuc(dna):
              counts = {'A': 0, 'C':0, 'G':0,'T':0}
              for nuc in dna:
                  counts[nuc] += 1
              return counts

                      DNAstring = "CCATTCCGTCCGGCCCTCAGCGTGAAACTCCACCTATACGCAAATATCCGCCGAATATAGCGATAAGACCAGTTCGGAGTAATCACTAGAAAAGGACGGCGTGTATCCTCGAGATGCAGATACCAGGTTCATCACTCTAGAAGACCCCATTCGTATGGCAGCCGTATACGCTGTTGCGTGATGCGCGTGCACAGATTTATCCGTACTTGCGGGATTGTCCCGGACCCCACTCCATCCATTCCTCGAAGCCAGGCAGACCTAGCGGTAGCGGGAGGATAGCCCCTAATGGCCCCTCTGTTTGCGCCCGATTAAGGGTCAGATCATTCTAGGTCAGAGAAAGTACATACCCTGTCCCTATACTTTCTACAGCCCCAGCCTGAAAGAATACCAGTAGACTTACTCGTATCCTTCCGAAGTTATGCTCACAGGGCAGGCTTGTGGGTACCACGAGCAGATTGTCGGCATATTATCACCCGCGCACTGATCCACCACGCACGCTCGTGTTGCGATCAATTGTTATAGTGGACCCATGTTGTATCCTAGCTCATCAATAACTTTGGTAAGCGCCAGAAGAGTGCATAGAACGCGGAATCACGATGAGATGCCAATCTACGATCTAGGTGGCACAATCTGGCCCTAACTGAACGACGCTGCTAACTTTATAAATTATCGTAGAGGGATAGACCTTACTAGTGGCTGGAAATACAATGTCTTTTGCGCTCTGTTCCACTGAATCCTAGTCCGGCCCACCACGCAATAATACGGACATAGGCGGAGACAGGGCGGTACCCGCGTTCGTACGTCTACCACGGTTGTAGGAGTCCTCGACTGGGATCCTCTGGGTCGGAGAGAAATGCTCTTCATCTTCCGTTCCATGTACATTGGCTTCATAGTATGCCCGCCGAGTTGCCCAGATCCGTCGTCTACCCTAACAACACAGCGGAGTTCGAATGTAGAAGTCCTTGA"

          result = countNuc(DNAstring)
          print(' '.join([str(val) for key, val in result.items()]))
          
          
          
          
          
          
EXERCISE 2 RNA


          def DNAtoRNA(seq):
                  return seq.replace("T", "U")
          print(DNAtoRNA("CAGTCACGAGAACACGTTCGCAGAACTAGGCGGTCGCTCAGCGCAACGACGCGACTTTCTCACAGGTTCCACTGTCCCCAGTGAGCTCCGAAAAAACCCTCCGCGTTGTCGATACGCCAACTGGTGCCGTCTCCGCCCATATCGGTATATTCTGTCCCGTACATGTAACATTTCTACACTCAACTCATGACTAGATCAGCGACTGACGCAGAGTTTACGTAGTAGACAGACCCAGCGGTAGGCTTTTCCTGTATCATGCATAACGTAGCCGCCATCCAGCGTAAAAATAAGCGGCAATTGCACAGACATCTCTGTTGATTATTTTGGAACGTTGCAATGGAACAGTAGCGCTTAATGTGAGTCACATGACTGTTCGTGACAGCACGGCGTGCCTGTTACTTATAGCAGGGCATCATACAGTCTTAGCGTCGTGTCCATCATTAACACCCTACTGCGTCTATGTCATTAGTAGCAGGATTCGTGCTAAGAGGGTCGTAAAGGCGTTTGTCTTTATAACCTTTTATCTTGTTGAAGGTACAGAATCGGGTAAGTCTATATCGTCACTCTTGGAATTAAAGCAGCATTACTGAGAGAAGCGGACATATTCCGGCACGAGTGCTCTAAACGCTTGGAAGAGTATATGGTCGGGGACACTCGACTATTCAACTGTTCCCCAACCAGACGAGTTGTTCTTTGCAATTCGATGGTCTCATGAAGTTGTCGAGGAAATTATTGTTCTCAACCAAGGGACACACAGCCGGGTCATAAGCCACTAAAAAACATTGACGATGTTCTCCTAAGTGATATCTCATTAAGCGGAAGTATGGGATATTGGCTCAGATATCTAGTGAACTACCTCATATAGCCTCCCCCATTCTGTTCGAACTAGAGTCGTAGCTGTACGAGTCTCACATAGCTATTAGCCGCATCGAGATCAGCTGGGAACCGCCGTCTTTATAACAGGTCTGG")
          
          
 EXERCISE 3 REVC
 
 
           DNA_reversecomplement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

          def reverse_complement(seq):
              return ''.join([DNA_reversecomplement[nuc] for nuc in seq])[::-1]
          print(reverse_complement("ATAGTCCCGCTACGGCGCCTGGGACCCATCGTAGACAGAGATGGAAAAGATGCCGGTCTACAATGAGCGGGTACCGTGGTGGCCGTTATCGCGCCCGTACAAGGACGGTCATTGCAGAGCGCCATACACTCCTGATGTGGCCGCGCCCATGATTGTATAGGTTTATTGCGGACTTGAATGGCGAAAAAGATGTCAACTTGGACGGCGATTGTCCGGCCTACAAGTATCACAGACGCCCCGGGGAGAGCTGTAAGCGAAGGGGATATCTCGCCTTCAATTCCTTCAGATCAAACCGTGTAGATAAAGATTGAAATTCTATGATGGGACATCCATGCCGAGACGGCAAGAGCATGGGTGATTCAGGCGCTCTTGCCCTAGTCGGAAGTTCTCCATCGTGTCAAAAGTTTAGGATCTGGCAGACCGTTCCTCTCCCAATCCCCCCTGGCTTATATTCTAGCGCACGGAACGGGTGTCCTTTATCATTTGCCGTCTACACGCCCACAGCTATCTGACCACCTAGTATCAAAAGTGCCATACTTGCCGTAAATCCTTAGCCCTCCGCTGACTTACATTCCATCTGGGCACGAATTCCCCACTATTCAAACACGAGGCGCGGTTTGCCTACCTCAACGCAGGTGGATTATGCACAATATTCAGTATACCTGAAAGGCCACGCGGAACTTCTGGGCATCACGCGTGTTCCTTACTGATTCTCGTGCCGGTGCCCTCCCCGACAGGGAAATCCTACACGGATTATGGCACGGTGGCCAATCCATCGAGGTTAGTCATATTGTGACAGGGGTTTGGTCGGGGATGTAAGGGAGTCAGAAGTTACGTCGCGTCCACTATACGGAGTTTTTAGCTGTGACCAGATTGAAACTAATGTAAC"))

