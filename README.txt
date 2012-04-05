For a complete documentation please go to http://archtk.github.com

Installation Requirements:

Python 2.6.x or 2.7.x
Numpy package http://numpy.scipy.org/
If you intend to use the vascular adaptation algorithm Scipy package is strongly recommended. (http://www.scipy.org)
A browser for visualizing simulation results (internet connection is not required). Currently supported browser are Firefox, Chrome, Safari, Opera and Internet Explorer

For additional features:
MEncoder library http://www.mplayerhq.hu/design7/dload.html (Only for velocity profile videos)
MatPlotLib library http://matplotlib.sourceforge.net/ (Post processing using png files instead of default browser visualization and for velocity profile videos.)
LXML package http://lxml.de/   (Only for xml validation feature).

How To:
Open a terminal and type:
./pyNS --help
For 

Open a terminal and type:
./pyNS
By default pyNS runs a vascular network which represents the arterial vasculature of a right arm.

Open results.html located in the Results folder for visualize results.

Guidelines for parameters.csv file, needed for arm templates simulations:
Date of birth (dob) has to be in dd/mm/yyyy
Date of surgery (dos) has to be in dd/mm/yyyy 
Gender 0 Female, 1 Male
Arm 0 Left, 1 Right
ftype 0 End-to-end radio-cephalic fistula, 1 End-to-side radio-cephalic fistula, 3 End-to-side brachio-cephalic fistula, 5 End-to-side brachio-basilic fistula
Diabetes 0 No, 1 Yes
Hypertension 0 No, 1 Yes
Height has to be in cm
Weight has to be in Kg
Systolic and diastolic pressures (says and diap) have to be in mmHg
Cardiac output is not mandatory but has to be in mL/min
Period is the cardiac period, measured in seconds
Brachial, radial and ulnar flows have to be in mL/min
Hematocrit (ht) has to be in %
Protein plasma concentration (cp) has to be in g/dL

If any parameter is not specified, pyNS assumes default values:
dos = 27/07/2010
dob = 27/07/1960
gender Male
Arm Right
Fistula Type End-to-side radio-cephalic
Height 175 cm
Weight 70 Kg
Systolic Pressure 110 mmHg
Diastolic Pressure 70 mmHg
Brachial Flow 130 mL/min
Radial Flow 20 mL/min
Ulnar Flow 30 mL/min
Period 1 second
Hematocrit 45 %
Protein Plasma Concentration 7 g/dL
Hypertension No
Diabetes No