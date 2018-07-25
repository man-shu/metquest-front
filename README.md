
# metquest-front
Web based user interface for MetQuest
## Installation Steps
* Clone the metquest-front directory, as follows:
```
git clone https://github.com/man-shu/metquest-front.git
```
* Install metquest via pip, as follows: 
```
sudo pip3 install metquest
```
or check [here](https://github.com/aarthi31/metquest-1) for other installation options.
* Install [flask](http://flask.pocoo.org/) and [flask-wtforms](http://flask-wtf.readthedocs.io/en/stable/) via pip, as follows:
```
sudo pip install flask
sudo pip install Flask-WTF
```
* Install [d3flux](https://github.com/pstjohn/d3flux), as follows:
  * Download and extract 'd3flux-master.zip' file from [here](https://github.com/pstjohn/d3flux).
  * Change to the extracted 'd3flux-master' directory and setup a development version of the same, as follows:
  
  ```
  cd d3flux-master
  python setup.py develop
  ```
  * Copy and replace the 'd3flux-master' directory with the one (having the same name) in the 'metquest-front' directory.
  
* Change to the 'metquest-front' directory and run either of the 'metquest-flask.py' or 'knockout.py', as follows:
```
cd metquest-front
python3 metquest-flask.py
```
or
```
cd metquest-front
python3 knockout.py
```
