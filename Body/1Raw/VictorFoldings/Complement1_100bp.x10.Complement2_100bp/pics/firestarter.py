import re
import string
import os


files=os.listdir("const")





for f in files: 

    # os.system('RNAfold -C </home/jester/VIC/sheeps/%s  -P DNA  > /home/jester/VIC/wolves/%s' % (f,f))
    # os.system('/home/jester/VIC/bin/./RNAfold --commands "/media/jester/Backup Plus/VIC/norcon/work.constraints" <"/media/jester/Backup Plus/VIC/norcon/const/%s" -P DNA  > "/media/jester/Backup Plus/VIC/norcon/result/%s"' % (f,f)) 
    
    os.system('/mnt/lustre/vshamanskiy/folding/RNAfold/./RNAfold --commands /mnt/lustre/vshamanskiy/folding/final/100com10com100/work.constraints </mnt/lustre/vshamanskiy/folding/final/100com10com100/const/%s -P DNA  > /mnt/lustre/vshamanskiy/folding/final/100com10com100/result/%s' % (f,f))
     
  