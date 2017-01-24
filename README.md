# larlitecv
Analysis processor framework for working with LArLite and LArCV data

## Getting Started

There are two examples available in the app folder: `Example` and `SelectionExample`.

You can also make a minimial application using the comamnd `gen_larlitecv_app.py`. For example,

    $ gen_larlitecv_app.py test


This will make a new folder, `app/test`, which provides

  * an empty class header, `test.h`
  * an empty class source, `test.cxx
  * an executable source file, `run_test.cxx`
  * a GNUmakfile
  * a LinkDef.h file (which helps ROOT build the symbols it needs to have your class be useable in the ROOT interpreter or callable in pyROOT)


