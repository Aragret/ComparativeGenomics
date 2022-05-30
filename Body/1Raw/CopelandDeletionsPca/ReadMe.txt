Dear Bill,

Thank you a lot for your reply and sorry for my long silence. 
I am about to resubmit my draft and try to compare numerically our results with your scores of deletions on the third principal component.
Is it possible for you to share with me your matrix, depicted in figure 4C (The 3rd PCA component weight vector for each element of the matrix)?


This is a quick and dirty compilation of all three component vectors (each 80x80 bins). # I guess he wanted to write 207 (see R code)
Format:
Component 0:
       Bin    Loading/Weight
       5592:5799,13878:14085 0.5945063738756831
The bin name is in the form of a:b,c:d, where
a = start of bin containing first deleted base pair
b = end of bin containing first deleted base pair
c = start of bin containing last deleted base pair
d = end of bin containing last deleted base pair
Note: when c>a, the deletion crosses the first reference base (position 1) # I guess he wanted to write "a>c"
