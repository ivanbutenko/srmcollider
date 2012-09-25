
# run in code folder
cp ../runcollider.py .
cp Release/ReleaseNotes . 
patch -p1 < Release/Release_patch1.patch
patch -p1 < Release/Release_patch2.patch
python setup.py sdist --manifest-only
python setup.py sdist --formats=gztar,zip
rm runcollider.py 
rm ReleaseNotes  

