============
Installation
============

For the most up-to-date features, and for the freshest, newest bugs, I
recommend using the mercurial repository.  See below.

You can do a full install, but that requires a few pre-requisites.  If
you do not have those pre-requisites you can install p4 as a pure
Python module with limited functionality, sometimes that might be all
you need.

P4 wants a fairly recent Python,
but I have not tried it out with Python 3 yet -- I doubt that it would
work.  You will want Python built with readline (as most are).  To tell whether you
have Python built with readline, call up interactive Python, type in
something like 2 + 2, hit <return>, note that the answer is 4, and
then hit <up arrow>.  If you get 2 + 2 back, you have a Python with
readline; if you get gobbledygoop, you don't.


Preparations for the full install
=================================

I've installed it on Linux and Mac OS X.  In either case, you need to
have some basic C-language programming tools, including a compiler,
libraries, headers, and so on.   

**Preparations for the full install on the Mac**


You need Xcode, available from the apple web site.

You need to install numpy, a python module, and the gsl library (Gnu
Scientific Library).  

In the past, I installed dependencies with MacPorts by
doing::

    sudo port install py27-numpy
    sudo port install gsl-devel
    sudo port install py27-sphinx
    sudo port install py27-tkinter

And that was mostly ok.  The only problem was that the 4 sphinx
executables (sphinx-build et al) had `-27` on the end of the name.
See the Sphinx Makefile for making the Sphinx docs, where it says how
to use executables with non-standard names.  Or just change the names.

Lately (on Maverics) I've had better luck building my Python 2.7 with Homebrew http://brew.sh/.
 
**Preparations for the full install on Linux**

If the installation process complains about lack of Python.h, then you
need what on Ubuntu would be called 'python-dev' or some such, which
is not included in the default Ubuntu install. 

I have recently installed p4
on a 64-bit Ubuntu 14.04, and had to::

    sudo apt-get install python-numpy
    sudo apt-get install libgsl0-dev
    sudo apt-get install python-dev

If you don't already have mercurial,
you will need it as well, and sphinx for the html documentation::

    sudo apt-get install mercurial
    sudo apt-get install python-sphinx

And if you want to use the GUI tree-drawing::

    sudo apt-get install python-tk


Get the mercurial version
=========================

The best way to get p4 is by using mercurial.  

1. If you don't already have it, install mercurial on your computer.
   On ubuntu it is::

    sudo apt-get install mercurial

   The command for mercurial is ``hg``

2. Choose a good place to put the source code, and go there.  A new
   directory will be made there in the next step.
 
3. Clone it.  The repository is at
   http://code.google.com/p/p4-phylogenetics/  There, the google
   folks say::

    Get a local copy of the p4-phylogenetics repository with this command:
    hg clone https://code.google.com/p/p4-phylogenetics/


  Doing that will make a directory called ``p4-phylogenetics`` filled
  with the files of the repository.  If you want the directory to be
  called something else, follow the ``hg clone ...`` command above
  with a directory name, and then your clone will be put in a
  directory with that name instead of ``p4-phylogenetics``. 

  Whatever you call it, if the ``hg clone`` command works that directory
  will be created and filled with the latest p4 source code.


Installing it in-place
======================

My fave way of using the mercurial version of p4 is to use it
in-place rather than installing it with ``setup.py``.  To use it in-place,
you need to 

1. Add the p4 hg directory, eg ``/usr/local/src/P4Hg`` to your ``PYTHONPATH``

2. Add the p4 hg bin directory, eg ``/usr/local/src/P4Hg/bin`` to your ``PATH``

3. Build the ``pf`` module, installing it in-place

For example if you install it in your home directory, to add the p4
hg directory to your ``PYTHONPATH``, you might add something like the
following line to your ``~/.profile`` or ``~/.bash_profile``::

  export PYTHONPATH=$HOME/src/P4Hg

(depending on where your P4 hg directory is, and what it is called), or
you can add ::

  export PYTHONPATH=$PYTHONPATH:$HOME/src/P4Hg

if you already have a ``PYTHONPATH`` defined.

The second thing you will want to do is to add the location of the p4
script to your ``PATH``.  Similar to adjusting the ``PYTHONPATH``
above, you can add a line like this to your  ``~/.profile`` or ``~/.bash_profile``::

  export PATH=$PATH:$HOME/src/P4Hg/bin

depending on where your P4 hg directory is, and what it is called.

To build the ``pf`` module, say::

   python setup.py build_ext -i

It might actually work.  If it doesn't, note the error messages that
flew by.  The earliest error message is usually a clue.

Installing the html docs
========================

You will need `sphinx <http://sphinx.pocoo.org>`_ .  On Ubuntu, its::

   sudo apt-get install python-sphinx

Then go to ``share/sphinxdoc`` in the p4 source, and do::

    make html

And then you can open ``_build/html/index.html`` with your browser.

Or, the docs are online, at `<http://p4.nhm.ac.uk>`_


Updating from hg
=================

The best part of installing it in-place is that it makes it easy to
update.  Generally all you need to do is to go to the p4 hg directory
and say::

  hg pull
  hg update

That pair of commands is usually
sufficient.  Occasionally there may have been changes to the
C-language code in the ``pf`` module.  If that is the case (would you
be able to see those files as they are updated?), and you use the
``pf`` module
then you would need to do::

 python setup.py build_ext -i

You would also need to do that when you install it in-place for the
first time, or if you make any changes to the C-language code
yourself.  If you are not sure it is needed, its ok to do it anyway.


Installing scqdist, the sub-cubic quartet distance module
=========================================================

See the directory Qdist in the source, with its own instructions.


To see if it works
==================

If, in your shell, you are still in the same directory that you built
it from, go to some other directory, or the following test will not work.

To see if you can load the package, start up python and then::

    import p4

To see if the p4 script works, say (perhaps from a new terminal) to
your shell (not in interactive python)::

    p4 --help

(Once it gets installed, if everything went perfectly and it still
does not work, try it in a new shell, or maybe even restart your
terminal program to refresh your PATH and PYTHONPATH.)

.. _completion_on_the_mac:

Completion in MacOS 10.5 and 10.6
=================================

P4 has a simple but useful completion module (I like it enough to use
it for all my python work) but file completion in the python that
comes with Mac OS after 10.5 is broken.  To fix it, you can either
install a better Python (my preferred option), or, to partially fix it, you can, in a file
'~/.p4/interactive' (that is a text file called 'interactive' that is
put in a directory called '.p4' in your home directory) put a line
that says 'var.readlineUsesEditline = True' (no quotes).  More info
about this is found in the file p4/Var.py.

The last time I tried that was on my old Snow Leopard box; it comes with Python 2.6.1, which is a little
oldish but not too bad, and it comes with numpy.  So to build it all I
needed was to add gsl-devel with MacPorts.  To use completion, I
needed to set var.readlineUsesEditline = True as described above.  It
then gave me method name completion, and doc strings, but no method
sigs (ie the stuff inside the parentheses, ie the method args).
However, the doc strings had the method sigs, so it was not too bad.

.. _completion_oddness:

Bash completion oddness
=======================

You may try completion from bash, but odd things happen.  For example,
you might want to read in a file myDataFile.nex, so you say::

    p4 myD<tab>

but then instead of completion, you get::

   p4 myDTraceback (most recent call last):
  File "/path/to/p4/bin/p4", line 68, in <module>
    func.readFile(f)
  File "/path/to/p4/modules/p4/func.py", line 356, in readFile
    raise Glitch, gm
  p4.Glitch.Glitch: 

  func.readFile(help)
      Can't open help.  Are you sure you have the right name?

This oddness is because you have bash_completion, and there is a
completion file for another p4 (from Perforce).  It would be found in
``/etc/bash_completion.d/`` on Ubuntu, or maybe ``/opt/local/etc/completion.d``
from MacPorts.  Well, assuming that you don't actually use that other
p4, you can remove that file, and that gets rid of the oddness.

.. Making an RPM
.. =============
.. I've barely tested this, but it worked for me, long ago. YMMV.
.. To make an rpm (both source and binary), say::
..     python setup.py bdist_rpm
.. To install the resulting binary rpm in the default location, say as
   root::
..     rpm -ivh p4-0.xx-1.ix86.rpm
.. If you didn't use an rpm to install your current python or gsl, so rpm
.. does not know that it exists, you might have to say as root::
..     rpm -ivh --nodeps p4-0.xx-1.ix86.rpm


Deinstallation
==============

.. If its an rpm, easy::
..   rpm -e p4

There is a func.uninstall() function, which may work.  You may need to
run it as root, or use sudo.

If that does not work, then recall that things get installed in 3
places.  Search out the Python package, the p4 script, and the
examples.



Installing p4 using setup.py
============================

This is the usual way that Python packages are installed, and is an
alternative to installing p4 in-place as described above.  It can be
done from the hg download.

If you are upgrading, you can un-install the previous version with the
p4 func.uninstall() function.  Depending on how it was installed, you
may need to be root or use sudo to do that.

Maybe you are starting with a downloaded svn version, or maybe you are
starting with the file p4-0.xx.tar.gz.  If the latter, unpack it in
your favourite source directory.  In the newly-created directory note
the file setup.py.  That file controls the build and installation.  It
installs 3 things:

    1.  **The p4 package.**          Goes where 3rd party packages go
                                Eg /usr/local/lib/python2.6/site-packages/

    2.  **The p4 script.**           Goes somewhere in your path
                                Eg /usr/local/bin

    3.  **The examples and info.**  Goes in a share/doc directory
                                Eg /usr/local/share/doc/

Simple install
--------------

First you can build it, without installing it, by saying::

    python setup.py build

(no need to be root, or use sudo, for the above step)

After building it, you then install it.  The default location for
installation is where python libraries are installed, and you as
JoeUser may not have file-writing permission to put files there, so
you may need to be root or use sudo for the next step.  Eg if you sudo
it, you can say::

    sudo python setup.py install


Installation variations
-----------------------

To get an option reminder, do::
    
    python setup.py install --help. 

To install it in your home directory, say::

    python setup.py install --home=~

If you install it in your home directory, 
    
- there is no need to be root or to use sudo

- if you do this you may need to setenv your PYTHONPATH to eg
        ~/lib/python.  Eg in your ~/.bash_profile you can put the
        line::

          export PYTHONPATH=$HOME/lib/python

- you may also need to set your PATH environment variable to
      ~/bin.  In many cases this will already be done, but if it is
      not, and you are using the bash shell, you can do something like::

          PATH=$PATH:$HOME/bin

      and then, after all your paths have been set, you should have a::

          export PATH



Where things go
---------------


The default installation location has a "root", which might be /usr or
/usr/local, or your home directory.

The default location for installation of the modules is something like::

    /usr/lib/python2.7/site-packages, or
    ~/lib/python

depending on the "root" of the installation, of course.

The default location for the script p4 is something like::

    /usr/bin

The default location for the examples is something like::

    /usr/share/doc/p4-0.xx/Examples

 
If you want to statically link your gsl libs
============================================

For those who may not want to do the usual dynamic linking of gsl
libs, it is possible to statically link the gsl libs to the pf.so
module when you build it.  See the ``setup.py``
file, and uncomment and adjust the ``extra_link_args`` line.



