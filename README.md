pyviral
=======

Python code to numerically integrate mathematical models spread of disease.

References:

Kermack, W. O. and McKendrick, A. G. "A Contribution to the
Mathematical Theory of Epidemics." Proc. Roy. Soc. Lond. A 115,
700-721, 1927.

https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology


I used a precursor to this code in the following blog posts:

http://valeriecoffman.com/blog/walking-dead-mathematics
https://datacommunitydc.squarespace.com/blog/2013/01/better-science-of-viral-marketing-part-2
https://datacommunitydc.squarespace.com/blog/2013/02/better-science-of-viral-marketing-part-3


Usage
-----

```
import pyviral
pyviral.run(1e6, 100, 2e-6)
```
