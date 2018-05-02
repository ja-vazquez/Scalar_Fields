

import os, sys

# python archivo.py    nombre_in.txt ymin ymax 
# python archivo_2.py  nombre_in y_name


cmd=""" 
python compute_contours.py Quintom_contour.txt -0.5 2
python plot_fgivenx.py Quintom_contour \\\\w_{DE} 
"""


os.system(cmd)

