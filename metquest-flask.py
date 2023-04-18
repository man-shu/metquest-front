from flask import Flask, flash, redirect, render_template, request, url_for
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
from datetime import datetime

ALLOWED_EXTENSIONS = {'xml'}

app = Flask(__name__)
app.config.from_object(__name__)
app.config['SECRET_KEY'] = '7d441f27d441f27567d441f2b6176a'
app.config['UPLOAD_PATH'] = os.path.join(os.getcwd(), datetime.now().strftime('FILES_%d-%m-%Y_%I:%M:%S'))
if not os.path.exists(os.path.join(app.config['UPLOAD_PATH'])):
    os.makedirs(os.path.join(app.config['UPLOAD_PATH']))

class MetQuest(Form):
    length = TextField('Cut-off:', validators=[validators.required()])
    source = TextField('Source Metabolites:', validators=[validators.required()])
    seeds = TextField('Seed Metabolites:', validators=[validators.required()])
    target = TextField('Target Metabolite:', validators=[validators.required()])

    def reset(self):
        blankData = MultiDict([ ('csrf', self.reset_csrf() ) ])
        self.process(blankData)

pathsHTML = []
pathsUnion = []
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
    mqForm = MetQuest(request.form)

    return render_template('metUpload.html', uploadForm = uploadForm, mqForm = mqForm )

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/uploader', methods = ['GET', 'POST'])
def uploader():
    uploadForm = UploadForm(request.form)
    global modids
    global fnames
    if request.method == 'POST' and 'file' in request.files:
        listof = request.files.getlist('file')
        if len(listof) == 0:
           flash('Error : No selected file')
           print(uploadForm.errors)
           return render_template('metUpload.html')
        if len(listof) > 0:
            mets_in_model = []
            j = 0
            lentrack = 0
            modids = []
            fnames = []
            for f in listof:
                if allowed_file(f.filename):
                    filename = secure_filename(f.filename)
                    fnames.append(filename)
                    f.save(os.path.join(app.config['UPLOAD_PATH'], filename))
                    if modl := cobra.io.read_sbml_model(
                        os.path.join(app.config['UPLOAD_PATH'], filename)
                    ):
                        modids.append(modl.id)
                        mets_in_model.extend(
                            {'label': "", 'value': "", 'category': modl.id}
                            for _ in range(len(modl.metabolites))
                        )
                        while j < (lentrack + len(modl.metabolites)):
                            for metab in modl.metabolites:
                                mets_in_model[j]['label'] = f"{metab.name} ( {metab.id} )"
                                mets_in_model[j]['value'] = f"{modl.id} {metab.id}"
                                j = j+1
                        lentrack = lentrack + len(modl.metabolites)
                    else:
                        flash(f'Error : Model {filename} not valid. Upload a .xml model file.')
                        print(uploadForm.errors)
                        return render_template('metUpload.html')
                else:
                    flash(f'Error : Model {f.filename} not valid. Upload a .xml model file.')
                    print(uploadForm.errors)
                    return render_template('metUpload.html')
            flash('Model Uploaded')
            return render_template('metUpload.html', mets_in_model = mets_in_model, fnames = fnames)
    else:
        flash('Error : No file selected')
        print(uploadForm.errors)
        return render_template('metUpload.html')


@app.route('/full1', methods = ['GET', 'POST'])
def full1():
    global pathsHTML
    if request.method == 'POST':
        ind = request.form["full"]
        ind = int(ind)
        pathway1 = pathsHTML[ind]
        return render_template('path1.html', pathway1 = pathway1)

@app.route('/full2', methods = ['GET', 'POST'])
def full2():
    global pathsUnion
    if request.method == 'POST':
        ind = request.form["unFull"]
        ind = int(ind)
        pathway2 = pathsUnion[ind]
        return render_template('path2.html', pathway2 = pathway2)


def create_pathway(namemap, G, r):
    form = MetQuest(request.form)

    cofacs = ["ppi_m","pi_m","adp_m","atp_m","nh4_c","h_p","h_c","atp_c","adp_c","nadph_c","nadp_c","nad_c","nadh_c","h2o_m","co2_c","h_m","ade_c","amp_c","h2o_c","pi_c","pi_p","ade_p"]
    global modids
    rxn = [namemap[r[a]] for a in range(len(r))]
    mets = []
    d = {"reactions":[{"id":"R1","name":"","metabolites":{},"lower_bound":0,"upper_bound":1000,"gene_reaction_rule":"","notes":{"map_info":{"reversibility":False,"hidden":False,"cofactors":{}}}}],"metabolites":[{"id":"A","name":"","compartment":"","notes":{"map_info":{"display_name":"A"}}}],"genes":[],"id":"simple_model","compartments":{},"notes":{"map_info":{}}}

    for _ in range(len(r)-1):
        d["reactions"].append({"id":"R1","name":"","metabolites":{},"lower_bound":0,"upper_bound":1000,"gene_reaction_rule":"","notes":{"map_info":{"flux":9,"reversibility":False,"hidden":False,"cofactors":{}}}})

    for a in range(len(rxn)):
        d["reactions"][a]["id"] = r[a].replace(" ", "")
        d["reactions"][a]["notes"]["map_info"]["display_name"] = rxn[a]
        preds = list(G.predecessors(r[a]))
        succs = list(G.successors(r[a]))
        for pred in preds:
            mets.append(pred)
            d["reactions"][a]["metabolites"][pred.replace(" ", "")] = -1
            if len(pred.split(" ", 1)) == 1:
                if pred.split(" ", 1)[0] in cofacs:
                    d["reactions"][a]["notes"]["map_info"]["cofactors"][pred.replace(" ", "")] = {}
            elif pred.split(" ", 1)[1] in cofacs:
                d["reactions"][a]["notes"]["map_info"]["cofactors"][pred.replace(" ", "")] = {}
        for succ in succs:
            mets.append(succ)
            d["reactions"][a]["metabolites"][succ.replace(" ", "")] = 1
            if len(succ.split(" ", 1)) == 1:
                if succ.split(" ", 1)[0] in cofacs:
                    d["reactions"][a]["notes"]["map_info"]["cofactors"][succ.replace(" ", "")] = {}
            elif succ.split(" ", 1)[1] in cofacs:
                d["reactions"][a]["notes"]["map_info"]["cofactors"][succ.replace(" ", "")] = {}

    mets = list(set(mets))

    for _ in range(len(mets)-1):
        d["metabolites"].append({"id":"A","name":"","compartment":"","notes":{"map_info":{"display_name":"A"}}})

    for a in range(len(mets)):
        d["metabolites"][a]["id"] = mets[a].replace(" ", "")
        if len(modids) == 1:
            d["metabolites"][a]["notes"]["map_info"]["display_name"] = mets[a].split(" ",1)[1]
        else:
            d["metabolites"][a]["notes"]["map_info"]["display_name"] = mets[a]
        if len(mets[a].split(" ",1)) == 1:
            if mets[a].split(" ",1)[0] in cofacs:
                d["metabolites"][a]["notes"]["map_info"]["display_name"] = mets[a].split(" ",1)[1]
                d["metabolites"][a]["notes"]["map_info"]["hidden"] =True
        elif mets[a].split(" ",1)[1] in cofacs:
            d["metabolites"][a]["notes"]["map_info"]["display_name"] = mets[a].split(" ",1)[1]
            d["metabolites"][a]["notes"]["map_info"]["hidden"] =True
    return d

def print_summary(scope, tar_met, pathways, cut_len, cyclic, namemap, src_met, seed_met, G, all_tar2src, diff_tar2src):
#    global pathway1
#    pathway1 = create_pathway(seed_met, cut_len, src_met, tar_met, 0, most_diff, namemap, G)
#    global pathway2
#    pathway2 = create_pathway(seed_met, cut_len, src_met, tar_met, 1, most_diff, namemap, G)
    global pathsHTML
    global pathsUnion
    pathsHTML = []
    labels = []
    pathsUnion = []
    labelsUnion = []
    labelsSummary = []
    len_scopes = []
    cut_branched_pathss = []
    all_cnts = []
    cut_cnts = []
    cut_cyclics = []
    minstepss = []
    for i in range(len(tar_met)):
        for j in range(len(src_met)):
            labelSummary = "%s %s to %s summary"%("_".join(modids), src_met[j], tar_met[i])
            labelsSummary.append(labelSummary)
            if len(diff_tar2src[tar_met[i]][src_met[j]]) == 0:
                if len(all_tar2src[tar_met[i]][src_met[j]]) == 0:
                    print('No pathway found for %s to %s'%(src_met[j],tar_met[i]))
                    pathsHTML.append(None)
                    label  = "%s %s to %s no pathway"%("_".join(modids), src_met[j], tar_met[i])
                    labels.append(label)

                    pathsUnion.append(None)
                    labelsUnion.append(label)

                    len_scopes.append(len(scope))
                    cut_branched_paths = len(all_tar2src[tar_met[i]][src_met[j]])
                    cut_branched_pathss.append(cut_branched_paths)
                    all_cnts.append(None)
                    cut_cnts.append(None)
                    cut_cyclics.append(None)
                    minstepss.append(None)
                else:
                    print('just one pathway found for %s to %s'%(src_met[j],tar_met[i]))
                    k = 0
                    r = all_tar2src[tar_met[i]][src_met[j]][k]
                    d = create_pathway(namemap, G, r)
                    jsonFileName = '%s_%s to %s_just one pathway.json'%("_".join(modids), src_met[j], tar_met[i])
                    with open(os.path.join(app.config['UPLOAD_PATH'], jsonFileName), 'w') as fp:
                        json.dump(d, fp, separators=(',',':'))
                    mod = cobra.io.json.load_json_model(os.path.join(app.config['UPLOAD_PATH'], jsonFileName))
                    page = flux_map(mod, display_name_format=lambda x: str(x.id), figsize=(1024,500),flux_dict={rxn.id: None for rxn in mod.reactions})
                    pathsHTML.append(page)
                    label  = jsonFileName
                    labels.append(label)
                    pathsUnion.append(page)
                    labelsUnion.append(label)

                    pred = G.predecessors
                    succ = G.successors
                    all_pathways_count = []
                    path_count = []
                    cyclic_pathway_count = []
                    all_reactions_involved = []
                    len_scope = len(scope)
                    len_scopes.append(len_scope)
                    cut_branched_paths = len(all_tar2src[tar_met[i]][src_met[j]])
                    cut_branched_pathss.append(cut_branched_paths)
                    for plen in pathways[tar_met[i]]:
                        all_pathways_count.append(len(pathways[tar_met[i]][plen]))
                        if plen <= int(cut_len):
                            path_count.append(len(pathways[tar_met[i]][plen]))
                        for p in pathways[tar_met[i]][plen]:
                            for reactions in p:
                                all_reactions_involved.append(reactions)
                    if tar_met[i] in cyclic:
                        for plen in cyclic[tar_met[i]]:
                            if plen <= int(cut_len):
                                cyclic_pathway_count.append(len(cyclic[tar_met[i]][plen]))
                    all_cnt = sum(all_pathways_count)
                    all_cnts.append(all_cnt)
                    cut_cnt = sum(path_count)
                    cut_cnts.append(cut_cnt)
                    cut_cyclic = sum(cyclic_pathway_count)
                    cut_cyclics.append(cut_cyclic)
                    minsteps = min(pathways[tar_met[i]])
                    minstepss.append(minsteps)
            else:
                wt = {}
                for a in range(len(all_tar2src[tar_met[i]][src_met[j]])):
                    for b in range(len(all_tar2src[tar_met[i]][src_met[j]][a])):
                        wt[all_tar2src[tar_met[i]][src_met[j]][a][b].replace(" ","")] = 0
                for a in range(len(all_tar2src[tar_met[i]][src_met[j]])):
                    for b in range(len(all_tar2src[tar_met[i]][src_met[j]][a])):
                        wt[all_tar2src[tar_met[i]][src_met[j]][a][b].replace(" ","")] += 1
                rUnion = set()
                for a in range(len(all_tar2src[tar_met[i]][src_met[j]])):
                    rUnion.update(all_tar2src[tar_met[i]][src_met[j]][a])
                rUnion = list(rUnion)
                dUnion = create_pathway(namemap, G, rUnion)
                jsonFileName = '%s_%s to %s_UNION.json'%("_".join(modids), src_met[j], tar_met[i])
                with open(os.path.join(app.config['UPLOAD_PATH'], jsonFileName), 'w') as fp:
                        json.dump(dUnion, fp, separators=(',',':'))
                modUnion = cobra.io.json.load_json_model(os.path.join(app.config['UPLOAD_PATH'], jsonFileName))
                solution = modUnion.optimize()
                a = 0
                for reacts in modUnion.reactions:
                    rs = str(reacts.id)
                    solution.fluxes[a] = wt[rs]
                    a = a+1
                pageUnion = flux_map(modUnion, display_name_format=lambda x: str(x.id), figsize=(1024,500),flux_dict= solution)
                pathsUnion.append(pageUnion)
                labelUnion  = jsonFileName
                labelsUnion.append(labelUnion)

                pred = G.predecessors
                succ = G.successors
                all_pathways_count = []
                path_count = []
                cyclic_pathway_count = []
                all_reactions_involved = []
                len_scope = len(scope)
                len_scopes.append(len_scope)
                cut_branched_paths = len(all_tar2src[tar_met[i]][src_met[j]])
                cut_branched_pathss.append(cut_branched_paths)
                for plen in pathways[tar_met[i]]:
                    all_pathways_count.append(len(pathways[tar_met[i]][plen]))
                    if plen <= int(cut_len):
                        path_count.append(len(pathways[tar_met[i]][plen]))
                    for p in pathways[tar_met[i]][plen]:
                        for reactions in p:
                            all_reactions_involved.append(reactions)
                if tar_met[i] in cyclic:
                    for plen in cyclic[tar_met[i]]:
                        if plen <= int(cut_len):
                            cyclic_pathway_count.append(len(cyclic[tar_met[i]][plen]))
                all_cnt = sum(all_pathways_count)
                all_cnts.append(all_cnt)
                cut_cnt = sum(path_count)
                cut_cnts.append(cut_cnt)
                cut_cyclic = sum(cyclic_pathway_count)
                cut_cyclics.append(cut_cyclic)
                minsteps = min(pathways[tar_met[i]])
                minstepss.append(minsteps)


                for k in range(2):
                    print('Most different pathway combination found for %s to %s'%(src_met[j],tar_met[i]))
                    r = diff_tar2src[tar_met[i]][src_met[j]][k]
                    d = create_pathway(namemap, G, r)
                    jsonFileName = '%s_%s to %s_%i.json'%("_".join(modids), src_met[j], tar_met[i], k)
                    with open(os.path.join(app.config['UPLOAD_PATH'], jsonFileName), 'w') as fp:
                        json.dump(d, fp, separators=(',',':'))
                    mod = cobra.io.json.load_json_model(os.path.join(app.config['UPLOAD_PATH'], jsonFileName))
                    page = flux_map(mod, display_name_format=lambda x: str(x.id), figsize=(1024,500),flux_dict={rxn.id: None for rxn in mod.reactions})
                    pathsHTML.append(page)
                    label  = jsonFileName
                    labels.append(label)
    return render_template('combi.html', pathsHTML = pathsHTML, labels = labels, pathsUnion = pathsUnion, labelsUnion = labelsUnion, labelsSummary = labelsSummary, len_scopes = len_scopes, cut_branched_pathss = cut_branched_pathss, all_cnts = all_cnts, cut_cnts = cut_cnts, cut_cyclics = cut_cyclics, minstepss = minstepss, cut_len = cut_len, tar_met = tar_met)

#    return render_template('summary.html',len_scope = len_scope, cut_branched_paths = cut_branched_paths, all_cnt = all_cnt, cut_cnt = cut_cnt, cut_cyclic = cut_cyclic, minsteps = minsteps, pathway1 = pathway1, pathway2 = pathway2, src_met = src_met, tar_met = tar_met, cut_len = cut_len )



@app.route("/test", methods=['GET', 'POST'])
def test():
    mqForm = MetQuest(request.form)
    if request.method == 'POST':
    #if request.form.post['action'] == 'make_paths':
       if mqForm.validate():
          global modids
          if len(modids) != 0 :
             cut_len = request.form['length']
             cut_len = int(cut_len)
             src_met = request.form['source']
             src_met = src_met.split(",")
             for i in range(len(src_met)):
                 if src_met[i] == "":
                    del src_met[i]
             seed_met = request.form['seeds']
             seed_met = seed_met.split(",")
             seed_met = set(seed_met)
             tar_met = request.form['target']
             tar_met = tar_met.split(",")
             for i in range(len(tar_met)):
                 if tar_met[i] == "":
                    del tar_met[i]
             G, namemap = mq.create_graph(os.path.join(app.config['UPLOAD_PATH']),len(modids))
             print('Graph made')
             pathways, cyclic, scope = mq.find_pathways(G,seed_met,cut_len)
             print('pathways found')
             diff_tar2src = {}
             all_tar2src = {}
             check1 = 0
             check2 = 0
             for i in range(len(tar_met)):
                 diff_tar2src[tar_met[i]] = {}
                 all_tar2src[tar_met[i]] = {}
                 for j in range(len(src_met)):
                     src = []
                     src.append(src_met[j])
                     most_diff, s2t = mq.find_pathways_starting_from_source(src, pathways, tar_met[i], cut_len, G)
                     if src_met[j] not in most_diff:
                        most_diff[src_met[j]] = {}
                     diff_tar2src[tar_met[i]][src_met[j]] = most_diff[src_met[j]]
                     all_tar2src[tar_met[i]][src_met[j]] = s2t
                     check1 = check1 + len(diff_tar2src[tar_met[i]][src_met[j]])
                     check2 = check2 + len(all_tar2src[tar_met[i]][src_met[j]])
             if check1 == 0 and check2 == 0:
                flash('Error : No pathways could be found. Consider changing the cut-off or the seed metabolite set.')
                print(mqForm.errors)
                return render_template('metUpload.html', mqForm = mqForm)
             else:
                print("creating path representations")
                return(print_summary(scope, tar_met, pathways, cut_len, cyclic, namemap, src_met, seed_met, G, all_tar2src, diff_tar2src))

          else:
             flash('Error : Model file not uploaded. Remember to upload after selecting the .xml model file.')
             print(mqForm.errors)
             return render_template('metUpload.html', mqForm = mqForm)

       else:
          flash('Error : All the form fields are required. ')
          print(mqForm.errors)
          return render_template('metUpload.html', mqForm = mqForm)

if __name__ == "__main__":
    url = 'http://127.0.0.1:5000'
    webbrowser.open_new(url)
    app.run()
