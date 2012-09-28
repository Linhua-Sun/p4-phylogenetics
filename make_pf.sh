# This would be useful for people who use the svn repository, and want
# to keep it up to date.  The idea is to use it where it is, rather
# than 'installing' it as would be usual for python modules.

# I also include this one-line script for anybody who might want to
# muck with p4 in the privacy of their own home directory, without
# actually installing the thing.  If this applies to you, don't
# install p4 via 'python setup.py install'.  Instead, simply unpack
# the source somewhere, and add the path to the unpacked source to
# your PYTHONPATH, eg

# export PYTHONPATH=$HOME/P4

# Then, after building pf (see below), you can import p4, and python
# knows where to look.

# If in addition you want to use the p4 script, you can add the p4/bin
# directory to your PATH, or alias it.

# But in order for the above to work you need to build the pf module,
# which is written in the c-language.  You can build it and move pf.so
# to the p4 directory using this --

rm -f p4/pf.so

python setup.py build_ext -i
#env ARCHFLAGS='-arch x86_64' python setup.py build_ext -i
#env ARCHFLAGS='-arch ppc -arch i386' python setup.py build_ext -i # This should not be needed, should be default.

# The c-language stuff does not change much over time, so this need
# not be done every update.  But occasionally something in the pf
# module changes, and so this build_ext will need to be done again.
