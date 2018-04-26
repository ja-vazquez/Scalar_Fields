

import os, sys

# python archivo.py    nombre_in.txt ymin ymax 
# python archivo_2.py  nombre_in y_name


cmd=""" 
python compute_contours.py pto_1.txt 0.0 1.8
python plot_fgivenx.py pto_1 \\\\rho_{DE} 
"""


os.system(cmd)

