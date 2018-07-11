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


ALLOWED_EXTENSIONS = set(['xml'])

app = Flask(__name__)
app.config.from_object(__name__)
app.config['SECRET_KEY'] = '7d441f27d441f27567d441f2b6176a'
app.config['UPLOAD_PATH'] = '/home/himanshu/Desktop/uploads/'

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


#@app.route('/full1', methods = ['GET', 'POST'])
#def full1():
#    global pathsHTML
#    if request.method == 'POST':
#        ind = request.form["full"]
#        ind = int(ind)
#        pathway1 = pathsHTML[ind]
#        return render_template('path1.html', pathway1 = pathway1)

#@app.route('/full2', methods = ['GET', 'POST'])
#def full2():
#    global pathsUnion
#    if request.method == 'POST':
#        ind = request.form["unFull"]
#        ind = int(ind)
#        pathway2 = pathsUnion[ind]
#        return render_template('path2.html', pathway2 = pathway2)


#def create_pathway(namemap, G, r):
#    form = MetQuest(request.form)

#    cofacs = ["ppi_m","pi_m","adp_m","atp_m","nh4_c","h_p","h_c","atp_c","adp_c","nadph_c","nadp_c","nad_c","nadh_c","h2o_m","co2_c","h_m","ade_c","amp_c","h2o_c","pi_c","pi_p","ade_p"]
#    global modids
#    rxn = []
#    for a in range(len(r)):
#        rxn.append(namemap[r[a]])
#    mets = []
#    d = {"reactions":[{"id":"R1","name":"","metabolites":{},"lower_bound":0,"upper_bound":1000,"gene_reaction_rule":"","notes":{"map_info":{"reversibility":False,"hidden":False,"cofactors":{}}}}],"metabolites":[{"id":"A","name":"","compartment":"","notes":{"map_info":{"display_name":"A"}}}],"genes":[],"id":"simple_model","compartments":{},"notes":{"map_info":{}}}

#    for a in range(len(r)-1):
#        d["reactions"].append({"id":"R1","name":"","metabolites":{},"lower_bound":0,"upper_bound":1000,"gene_reaction_rule":"","notes":{"map_info":{"flux":9,"reversibility":False,"hidden":False,"cofactors":{}}}})

#    for a in range(len(rxn)):
#        d["reactions"][a]["id"] = r[a].replace(" ", "")
#        d["reactions"][a]["notes"]["map_info"]["display_name"] = rxn[a]
#        preds = list(G.predecessors(r[a]))
#        succs = list(G.successors(r[a]))
#        for b in range(len(preds)):
#            mets.append(preds[b])
#            d["reactions"][a]["metabolites"][preds[b].replace(" ", "")] =-1
#            if len(preds[b].split(" ",1)) == 1:
#                if preds[b].split(" ",1)[0] in cofacs:
#                    d["reactions"][a]["notes"]["map_info"]["cofactors"][preds[b].replace(" ", "")] ={}
#            else :
#                if preds[b].split(" ",1)[1] in cofacs:
#                    d["reactions"][a]["notes"]["map_info"]["cofactors"][preds[b].replace(" ", "")] ={}
#        for c in range(len(succs)):
#            mets.append(succs[c])
#            d["reactions"][a]["metabolites"][succs[c].replace(" ", "")] =1
#            if len(succs[c].split(" ",1)) == 1:
#                if succs[c].split(" ",1)[0] in cofacs:
#                    d["reactions"][a]["notes"]["map_info"]["cofactors"][succs[c].replace(" ", "")] ={}
#            else :
#                if succs[c].split(" ",1)[1] in cofacs:
#                    d["reactions"][a]["notes"]["map_info"]["cofactors"][succs[c].replace(" ", "")] ={}

#    mets = list(set(mets))

#    for a in range(len(mets)-1):
#        d["metabolites"].append({"id":"A","name":"","compartment":"","notes":{"map_info":{"display_name":"A"}}})

#    for a in range(len(mets)):
#        d["metabolites"][a]["id"] = mets[a].replace(" ", "")
#        if len(modids) == 1:
#            d["metabolites"][a]["notes"]["map_info"]["display_name"] = mets[a].split(" ",1)[1]
#        else:
#            d["metabolites"][a]["notes"]["map_info"]["display_name"] = mets[a]
#        if len(mets[a].split(" ",1)) == 1:
#            if mets[a].split(" ",1)[0] in cofacs:
#                d["metabolites"][a]["notes"]["map_info"]["display_name"] = mets[a].split(" ",1)[1]
#                d["metabolites"][a]["notes"]["map_info"]["hidden"] =True
#        else:
#            if mets[a].split(" ",1)[1] in cofacs:
#                d["metabolites"][a]["notes"]["map_info"]["display_name"] = mets[a].split(" ",1)[1]
#                d["metabolites"][a]["notes"]["map_info"]["hidden"] =True
#    return d

#def print_summary(scope, tar_met, pathways, cut_len, cyclic, namemap, src_met, seed_met, G, all_tar2src, diff_tar2src):
##    global pathway1
##    pathway1 = create_pathway(seed_met, cut_len, src_met, tar_met, 0, most_diff, namemap, G)
##    global pathway2
##    pathway2 = create_pathway(seed_met, cut_len, src_met, tar_met, 1, most_diff, namemap, G)
#    global pathsHTML
#    global pathsUnion
#    pathsHTML = []
#    labels = []
#    pathsUnion = []
#    labelsUnion = []
#    labelsSummary = []
#    len_scopes = []
#    cut_branched_pathss = []
#    all_cnts = []
#    cut_cnts = []
#    cut_cyclics = []
#    minstepss = []
#    for i in range(len(tar_met)):
#        for j in range(len(src_met)):
#            labelSummary = "%s %s to %s summary"%("_".join(modids), src_met[j], tar_met[i])
#            labelsSummary.append(labelSummary)
#            if len(diff_tar2src[tar_met[i]][src_met[j]]) == 0:
#                if len(all_tar2src[tar_met[i]][src_met[j]]) == 0:
#                    print('No pathway found for %s to %s'%(src_met[j],tar_met[i]))
#                    pathsHTML.append(None)
#                    label  = "%s %s to %s no pathway"%("_".join(modids), src_met[j], tar_met[i])
#                    labels.append(label)

#                    pathsUnion.append(None)
#                    labelsUnion.append(label)

#                    len_scopes.append(len(scope))
#                    cut_branched_paths = len(all_tar2src[tar_met[i]][src_met[j]])
#                    cut_branched_pathss.append(cut_branched_paths)
#                    all_cnts.append(None)
#                    cut_cnts.append(None)
#                    cut_cyclics.append(None)
#                    minstepss.append(None)
#                else:
#                    print('just one pathway found for %s to %s'%(src_met[j],tar_met[i]))
#                    k = 0
#                    r = all_tar2src[tar_met[i]][src_met[j]][k]
#                    d = create_pathway(namemap, G, r)
#                    with open('/home/himanshu/Desktop/%s_%s to %s_just one pathway.json'%("_".join(modids), src_met[j], tar_met[i]), 'w') as fp:
#                        json.dump(d, fp, separators=(',',':'))
#                    mod = cobra.io.json.load_json_model('/home/himanshu/Desktop/%s_%s to %s_just one pathway.json'%("_".join(modids), src_met[j], tar_met[i]))
#                    page = flux_map(mod, display_name_format=lambda x: str(x.id), figsize=(1024,500),flux_dict={rxn.id: None for rxn in mod.reactions})
#                    pathsHTML.append(page)
#                    label  = "%s %s to %s just one pathway"%("_".join(modids), src_met[j], tar_met[i])
#                    labels.append(label)
#                    pathsUnion.append(page)
#                    labelsUnion.append(label)

#                    pred = G.predecessors
#                    succ = G.successors
#                    all_pathways_count = []
#                    path_count = []
#                    cyclic_pathway_count = []
#                    all_reactions_involved = []
#                    len_scope = len(scope)
#                    len_scopes.append(len_scope)
#                    cut_branched_paths = len(all_tar2src[tar_met[i]][src_met[j]])
#                    cut_branched_pathss.append(cut_branched_paths)
#                    for plen in pathways[tar_met[i]]:
#                        all_pathways_count.append(len(pathways[tar_met[i]][plen]))
#                        if plen <= int(cut_len):
#                            path_count.append(len(pathways[tar_met[i]][plen]))
#                        for p in pathways[tar_met[i]][plen]:
#                            for reactions in p:
#                                all_reactions_involved.append(reactions)
#                    if tar_met[i] in cyclic:
#                        for plen in cyclic[tar_met[i]]:
#                            if plen <= int(cut_len):
#                                cyclic_pathway_count.append(len(cyclic[tar_met[i]][plen]))
#                    all_cnt = sum(all_pathways_count)
#                    all_cnts.append(all_cnt)
#                    cut_cnt = sum(path_count)
#                    cut_cnts.append(cut_cnt)
#                    cut_cyclic = sum(cyclic_pathway_count)
#                    cut_cyclics.append(cut_cyclic)
#                    minsteps = min(pathways[tar_met[i]])
#                    minstepss.append(minsteps)
#            else:
#                wt = {}
#                for a in range(len(all_tar2src[tar_met[i]][src_met[j]])):
#                    for b in range(len(all_tar2src[tar_met[i]][src_met[j]][a])):
#                        wt[all_tar2src[tar_met[i]][src_met[j]][a][b].replace(" ","")] = 0
#                for a in range(len(all_tar2src[tar_met[i]][src_met[j]])):
#                    for b in range(len(all_tar2src[tar_met[i]][src_met[j]][a])):
#                        wt[all_tar2src[tar_met[i]][src_met[j]][a][b].replace(" ","")] += 1
#                rUnion = set()
#                for a in range(len(all_tar2src[tar_met[i]][src_met[j]])):
#                    rUnion.update(all_tar2src[tar_met[i]][src_met[j]][a])
#                rUnion = list(rUnion)
#                dUnion = create_pathway(namemap, G, rUnion)
#                with open('/home/himanshu/Desktop/%s_%s to %s_UNION.json'%("_".join(modids), src_met[j], tar_met[i]), 'w') as fp:
#                        json.dump(dUnion, fp, separators=(',',':'))
#                modUnion = cobra.io.json.load_json_model('/home/himanshu/Desktop/%s_%s to %s_UNION.json'%("_".join(modids), src_met[j], tar_met[i]))
#                solution = modUnion.optimize()
#                a = 0
#                for reacts in modUnion.reactions:
#                    rs = str(reacts.id)
#                    solution.fluxes[a] = wt[rs]
#                    a = a+1
#                pageUnion = flux_map(modUnion, display_name_format=lambda x: str(x.id), figsize=(1024,500),flux_dict= solution)
#                pathsUnion.append(pageUnion)
#                labelUnion  = "%s %s to %s UNION"%("_".join(modids), src_met[j], tar_met[i])
#                labelsUnion.append(labelUnion)

#                pred = G.predecessors
#                succ = G.successors
#                all_pathways_count = []
#                path_count = []
#                cyclic_pathway_count = []
#                all_reactions_involved = []
#                len_scope = len(scope)
#                len_scopes.append(len_scope)
#                cut_branched_paths = len(all_tar2src[tar_met[i]][src_met[j]])
#                cut_branched_pathss.append(cut_branched_paths)
#                for plen in pathways[tar_met[i]]:
#                    all_pathways_count.append(len(pathways[tar_met[i]][plen]))
#                    if plen <= int(cut_len):
#                        path_count.append(len(pathways[tar_met[i]][plen]))
#                    for p in pathways[tar_met[i]][plen]:
#                        for reactions in p:
#                            all_reactions_involved.append(reactions)
#                if tar_met[i] in cyclic:
#                    for plen in cyclic[tar_met[i]]:
#                        if plen <= int(cut_len):
#                            cyclic_pathway_count.append(len(cyclic[tar_met[i]][plen]))
#                all_cnt = sum(all_pathways_count)
#                all_cnts.append(all_cnt)
#                cut_cnt = sum(path_count)
#                cut_cnts.append(cut_cnt)
#                cut_cyclic = sum(cyclic_pathway_count)
#                cut_cyclics.append(cut_cyclic)
#                minsteps = min(pathways[tar_met[i]])
#                minstepss.append(minsteps)


#                for k in range(2):
#                    print('Most different pathway combination found for %s to %s'%(src_met[j],tar_met[i]))
#                    r = diff_tar2src[tar_met[i]][src_met[j]][k]
#                    d = create_pathway(namemap, G, r)
#                    with open('/home/himanshu/Desktop/%s_%s to %s_%i.json'%("_".join(modids), src_met[j], tar_met[i], k), 'w') as fp:
#                        json.dump(d, fp, separators=(',',':'))
#                    mod = cobra.io.json.load_json_model('/home/himanshu/Desktop/%s_%s to %s_%i.json'%("_".join(modids), src_met[j], tar_met[i], k))
#                    page = flux_map(mod, display_name_format=lambda x: str(x.id), figsize=(1024,500),flux_dict={rxn.id: None for rxn in mod.reactions})
#                    pathsHTML.append(page)
#                    label  = "%s %s to %s %i"%("_".join(modids), src_met[j], tar_met[i], k)
#                    labels.append(label)
#    return render_template('combi.html', pathsHTML = pathsHTML, labels = labels, pathsUnion = pathsUnion, labelsUnion = labelsUnion, labelsSummary = labelsSummary, len_scopes = len_scopes, cut_branched_pathss = cut_branched_pathss, all_cnts = all_cnts, cut_cnts = cut_cnts, cut_cyclics = cut_cyclics, minstepss = minstepss, cut_len = cut_len, tar_met = tar_met)

##    return render_template('summary.html',len_scope = len_scope, cut_branched_paths = cut_branched_paths, all_cnt = all_cnt, cut_cnt = cut_cnt, cut_cyclic = cut_cyclic, minsteps = minsteps, pathway1 = pathway1, pathway2 = pathway2, src_met = src_met, tar_met = tar_met, cut_len = cut_len )



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
             G, namemap = mq.create_graph("/home/himanshu/Desktop/uploads",len(modids))
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
    app.run(use_reloader=True, debug=True, port = 5002)
