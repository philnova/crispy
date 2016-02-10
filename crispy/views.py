#from flask import Flask, render_template, url_for, request, redirect, flash, jsonify
from flask import Flask, render_template, url_for, request, redirect, flash, jsonify
from crispy import app

from wtforms import Form, TextAreaField, RadioField, DateField, BooleanField, StringField, validators, FloatField, DecimalField

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from crispy.models import Base, Sequence

import datetime

import RawSequence

import jinja2
environment = jinja2.Environment() # you can define characteristics here, like telling it to load templates, etc


###########
# define forms
############

class BaseSearchForm(Form):
	species = RadioField(choices=[('GCA_000001405.15_GRCh38_no_alt_analysis_set', 'Human (GRCh38)'),('mm8','Mouse (UCSC mm8)'),('d_melanogaster_fb5_22','D. melanogaster (Flybase r5.22)'),('c_elegans_ws200','C. Elegans (Wormbase WS200)'),('s_cerevisiae','S. Cerevisiae (CYGD)'),('e_coli','E. Coli (NCBI st.536)')], default='GCA_000001405.15_GRCh38_no_alt_analysis_set')

class RawSequenceForm(BaseSearchForm):
	sequence = TextAreaField('Sequence', [validators.Required()], id="sequenceSearchField")

class CoordinateForm(BaseSearchForm):
	coordinate = StringField('Coordinate', [validators.Required()])

class GeneNameForm(BaseSearchForm):
	name = StringField('Gene Name', [validators.Required()])


#jinja2 tests
def high_quality(score):
	return score >= .75

def medium_quality(score):
	return .5 <= score < .75

def low_quality(score):
	return score < .5

environment.tests["high_quality"] = high_quality
environment.tests["medium_quality"] = medium_quality
environment.tests["low_quality"] = low_quality

#app = Flask(__name__)

engine = create_engine('sqlite:///guidesequences.db')
Base.metadata.bind = engine

DBSession = sessionmaker(bind=engine)
session = DBSession()

@app.route('/')
def mainPage():
	seq_form = RawSequenceForm(request.form)
	coord_form = CoordinateForm(request.form)
	name_form = GeneNameForm(request.form)
	return render_template('home.html', seq_form=RawSequenceForm(), coord_form=CoordinateForm(), name_form=GeneNameForm())

@app.route('/seqSearch', methods=['GET','POST'])
def seqSearch():
	seq_form = RawSequenceForm(request.form)
	if request.method == 'POST' and seq_form.validate():
		species = seq_form.species.data

		sequence_string = seq_form.sequence.data
		sequence_string.strip()
		sequence_string.replace(" ","")
		sequence_string = sequence_string.encode('utf-8')
		print sequence_string

		seq = RawSequence.RawSequence(sequence_string, genome=species)
		print seq

		results = seq.score_all()
		guides, offtargets, scores = [], [], []
		
		for res in results:
			guides.append(res[0])
			offtargets.append(res[1])
			scores.append(res[2])

		goodGuides, mediumGuides, badGuides = False, False, False #flags to notify page whether any such guides exist
		for score in scores:
			if goodGuides and mediumGuides and badGuides:
				break
			else:
				if score < 0.5:
					badGuides = True
				elif score < 0.75:
					mediumGuides = True
				else:
					goodGuides = True

		return render_template('seqResults.html', sequence=sequence_string, species=species, guides=guides,offtargets=offtargets, scores=scores, goodGuides = goodGuides, mediumGuides = mediumGuides, badGuides=badGuides)



@app.route('/coordSearch', methods=['GET','POST'])
def coordSearch():
	coord_form = CoordinateForm(request.form)
	if request.method == 'POST' and coord_form.validate():
		species = coord_form.species.data
		coordinate = coord_form.coordinate.data
		return redirect(url_for('coordResults', coordinate=coordinate, species=species))

@app.route('/coordResults/<string:species>/<string:coordinate>')
def coordResults(coordinate, species):
	return render_template('coordResults.html',coordinate=coordinate, species=species)

@app.route('/nameSearch', methods=['GET','POST'])
def nameSearch():
	name_form = GeneNameForm(request.form)
	if request.method == 'POST' and name_form.validate():
		species = name_form.species.data
		name = name_form.name.data
		return redirect(url_for('nameResults', name=name, species=species))

@app.route('/nameResults/<string:species>/<string:name>')
def nameResults(name, species):
	return render_template('nameResults.html',name=name, species=species)

if __name__ == '__main__':
    app.secret_key = 'super_secret_key' # gives us access to the session to flash messages
    app.debug = True
    app.run(host='0.0.0.0', port=5000)