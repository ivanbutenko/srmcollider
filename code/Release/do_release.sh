
# run in code folder
cp ../runcollider.py
python setup.py sdist --manifest-only
python setup.py sdist --formats=gztar,zip
patch -p1 < Release/Release_patch1.patch

