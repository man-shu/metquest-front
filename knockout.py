from flask import Flask, flash, redirect, render_template, request, url_for, send_file
from flask_wtf import Form, FlaskForm
from wtforms import Form, TextField, TextAreaField, validators, StringField, SubmitField, SelectField, IntegerField
from flask_wtf.file import FileField, FileAllowed, FileRequired
from werkzeug import secure_filename

import os
import metquest as mq
import json
import cobra
import webbrowser
from cobra.io import load_json_model
from cobra.core import Metabolite, Reaction, Model
from d3flux import flux_map
import matplotlib.pyplot as plot
from collections import Counter
import io
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from datetime import datetime

ALLOWED_EXTENSIONS = set(['xml'])

app = Flask(__name__)
app.config.from_object(__name__)
app.config['SECRET_KEY'] = '7d441f27d441f27567d441f2b6176a'
app.config['UPLOAD_PATH'] = os.path.join(os.getcwd(), datetime.now().strftime('KOFILES_%d-%m-%Y_%I:%M:%S'))
if not os.path.exists(os.path.join(app.config['UPLOAD_PATH'])):
    os.makedirs(os.path.join(app.config['UPLOAD_PATH']))
    print('file made')

class Knock(Form):
    length = TextField('Cut-off:', validators=[validators.required()])
    knock_out = TextField('Reaction to be knocked-out:', validators=[validators.required()])
    seeds = TextField('Seed Metabolites:', validators=[validators.required()])
    target = TextField('Target Metabolite:', validators=[validators.required()])

    def reset(self):
        blankData = MultiDict([ ('csrf', self.reset_csrf() ) ])
        self.process(blankData)

pathsHTML = []
cobra_mods = []
modids = []
fnames = []

class UploadForm(FlaskForm):
    upload = FileField(validators=[FileRequired()])

@app.route('/')
def inputs():
    global modids
    global fnames
    modids = []
    fnames = []
    uploadForm = UploadForm(request.form)
    knockForm = Knock(request.form)

    return render_template('knockUpload.html', uploadForm = uploadForm, knockForm = knockForm )

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/uploader', methods = ['GET', 'POST'])
def uploader():
   uploadForm = UploadForm(request.form)
   global modids
   global fnames
   global cobra_mods
   cobra_mods = []
   modids = []
   fnames = []
   if request.method == 'POST' and 'file' in request.files:
      listof = request.files.getlist('file')
      if len(listof) == 0:
         flash('Error : No selected file')
         print(uploadForm.errors)
         return render_template('knockUpload.html')
      if len(listof) > 0:
         rxn_in_model = []
         j = 0
         lentrack = 0
         for f in listof:
             if allowed_file(f.filename):
                filename = secure_filename(f.filename)
                fnames.append(filename)
                f.save(os.path.join(app.config['UPLOAD_PATH'], filename))
                modl=cobra.io.read_sbml_model(os.path.join(app.config['UPLOAD_PATH'], filename))
                cobra_mods.append(modl)
                if modl:
                   modids.append(modl.id)
                   for i in range(len(modl.reactions)):
                       rxn_in_model.append({'label':"", 'value':"", 'category': modl.id})
                   while j < (lentrack + len(modl.reactions)):
                         for reacts in modl.reactions:
                             rxn_in_model[j]['label'] = str(reacts.name) + " ( " + str(reacts.id) + " )"
                             rxn_in_model[j]['value'] = str(modl.id) + " " + str(reacts.id)
                             j = j+1
                   lentrack = lentrack + len(modl.reactions)
                else:
                   flash('Error : Model %s not valid. Upload a .xml model file.'%(filename))
                   print(uploadForm.errors)
                   return render_template('knockUpload.html')
             else:
                flash('Error : Model %s not valid. Upload a .xml model file.'%(f.filename))
                print(uploadForm.errors)
                return render_template('knockUpload.html')
         flash('Model Uploaded')
         return render_template('knockUpload.html', rxn_in_model = rxn_in_model, fnames = fnames)
   else:
      flash('Error : No file selected')
      print(uploadForm.errors)
      return render_template('knockUpload.html')

@app.route("/test", methods=['GET', 'POST'])
def test():
    knockForm = Knock(request.form)
    global modids
    if request.method == 'POST':
    #if request.form.post['action'] == 'make_paths':
       if knockForm.validate():
          global modids
          if len(modids) != 0 :
             cut_len = request.form['length']
             cut_len = int(cut_len)
             knockout = request.form['knock_out']
             knockout = list(set(knockout.split(",")))
             for i in knockout:
                 if i == "":
                    knockout.remove(i)
             seed_met = request.form['seeds']
             seed_met = seed_met.split(",")
             seed_met = set(seed_met)
             tar_met = request.form['target']
             tar_met = tar_met.split(",")
             for i in tar_met:
                 if i == "":
                    tar_met.remove(i)
             G, namemap = mq.create_graph(os.path.join(app.config['UPLOAD_PATH']),len(modids))
             print('Graph made')
             pathways, cyclic, scope = mq.find_pathways(G,seed_met,cut_len)
             print('pathways found')
             all_reactions = {}
             freq = {}
             check = 0
             for i in range(len(tar_met)):
                 all_reactions[tar_met[i]] = []
                 pred = G.predecessors
                 succ = G.successors
                 if tar_met[i] in pathways:
                    all_reactions_involved = []
                    for plen in pathways[tar_met[i]]:
                        for paths in pathways[tar_met[i]][plen]:
                            for reactions in paths:
                                all_reactions_involved.append(namemap[reactions])
                    all_reactions[tar_met[i]] = all_reactions_involved
                 freq[tar_met[i]] = dict(Counter(all_reactions[tar_met[i]]))
                 for keys,values in freq[tar_met[i]].items():
                    freq[tar_met[i]][keys] = values/len(all_reactions[tar_met[i]])
                 check = check + len(all_reactions[tar_met[i]])
             if check == 0 :
                flash('Error : No pathways could be found. Consider changing the cut-off or the seed metabolite set.')
                print(knockForm.errors)
                return render_template('knockUpload.html', knockForm = knockForm)
             else:
                print("knocking out reaction(s)")
                mods_to_knock_from = {}
                for i in knockout:
                    mods_to_knock_from[i.split(" ")[0]] = []
                for i in knockout:
                    mods_to_knock_from[i.split(" ")[0]].append(i.split(" ")[1])
                for keys, values in mods_to_knock_from.items():
                    for j in cobra_mods:
                        if keys == j.id:
                            j.remove_reactions(values, True)
                            Jid = j.id + "-" + "-".join(values)
                            modids = []
                            modids.append(Jid)
                fold_name = "_".join(modids)
                if not os.path.exists(os.path.join(app.config['UPLOAD_PATH'], fold_name)):
                    os.makedirs(os.path.join(app.config['UPLOAD_PATH'], fold_name))
                for i in cobra_mods:
                    cobra.io.write_sbml_model(i, os.path.join(app.config['UPLOAD_PATH'], fold_name + "/" + i.id + ".xml"))
                G_k, namemap_k = mq.create_graph(os.path.join(app.config['UPLOAD_PATH'], fold_name), len(modids))
                pathways_k, cyclic_k, scope_k = mq.find_pathways(G_k,seed_met,cut_len)
                all_reactions_k = {}
                freq_k = {}
                final_freq_k ={}
                check1 = 0
                for i in range(len(tar_met)):
                    all_reactions_k[tar_met[i]] = []
                    if tar_met[i] in pathways_k:
                        all_reactions_involved_k = []
                        for plen in pathways_k[tar_met[i]]:
                            for paths in pathways_k[tar_met[i]][plen]:
                                for reactions in paths:
                                    all_reactions_involved_k.append(namemap_k[reactions])
                        all_reactions_k[tar_met[i]] = all_reactions_involved_k
                    freq_k[tar_met[i]] = dict(Counter(all_reactions_k[tar_met[i]]))
                    for keys,values in freq_k[tar_met[i]].items():
                        freq_k[tar_met[i]][keys] = values/len(all_reactions_k[tar_met[i]])
                    final_freq_k[tar_met[i]] = {}
                    y1 = []
                    x = list(freq[tar_met[i]].keys())
                    y2 = []
                    for keys, values in freq[tar_met[i]].items():
                        final_freq_k[tar_met[i]][keys] = 0
                    for keys, values in freq_k[tar_met[i]].items():
                        final_freq_k[tar_met[i]][keys] = values
                    for react in x:
                        y1.append(freq[tar_met[i]][react])
                        y2.append(final_freq_k[tar_met[i]][react])
                    import numpy as np
                    x_nums = np.arange(1,len(freq[tar_met[i]])+1)
                    fig, ax = plot.subplots(figsize=(15,50))
                    p1 = ax.barh(x_nums + 0.25, y1, 0.25, align ='center', color="b", label='Before Knockout')
                    p2 = ax.barh(x_nums, y2, 0.25, align ='center',color= "r", label='After Knockout')
                    ax.set(yticks = x_nums + 0.25/2, yticklabels=x)
                    ax.legend()
                    ax.autoscale_view()
                    img = io.BytesIO()
                    plot.savefig(img)
                    img.seek(0)
                    check1 = check1 + len(all_reactions_k[tar_met[i]])
                return send_file(img, mimetype='image/png')

                
          else:
             flash('Error : Model file not uploaded. Remember to upload after selecting the .xml model file.')
             print(knockForm.errors)
             return render_template('knockUpload.html', knockForm = knockForm)

       else:
          flash('Error : All the form fields are required. ')
          print(knockForm.errors)
          return render_template('knockUpload.html', knockForm = knockForm)

if __name__ == "__main__":
    url = 'http://127.0.0.1:5002'
    webbrowser.open_new(url)
    app.run()
