Example code using Jupyter notebooks (formerly known as IPython notebooks)

To use this, first make sure your environment is setup properly:

* setup larlite (you should see the environemnt variable LARLITE_BASEDIR defined)
* setup larcv (you should see LARCV_BASEDIR defined)
* seutp larlitecv (likewise LARLITECV_BASEDIR)

Then run

    jupyter notebook


A web browser should appear. Go to example.ipynb. You can then follow along interactively.

# making your own notebook

If you want to make your own example, please don't edit the one saved in the repository.

Best thing to do is

* make a new branch of the repository: git checkout -b develop my_total_awesome_branch
* make a new folder in the app directory
* add that folder to the app directory
* and make a copy of example.ipynb

Now you can save your notebook.  
Also, you can check in your branch if you want (so others can look at it). 
This is one way to either share results, provide examples of how your code works, or recreate a bug for someone to help fix it.
