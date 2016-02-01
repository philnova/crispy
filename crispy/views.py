#from flask import Flask, render_template, url_for, request, redirect, flash, jsonify
from flask import Flask, render_template, url_for, request, redirect, flash, jsonify
from crispy import app

from wtforms import Form, TextAreaField, RadioField, DateField, BooleanField, StringField, validators, FloatField, DecimalField

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from crispy.models import Base, Sequence

import datetime

###########
# define forms
############

class RawSequenceForm(Form):
	sequence = TextAreaField('Sequence', [validators.Required()])

class CoordinateForm(Form):
	coordinate = StringField('Coordinate', [validators.Required()])

class GeneNameForm(Form):
	name = StringField('Gene Name', [validators.Required()])


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
		sequence = seq_form.sequence.data
		return redirect(url_for('seqResults', sequence=sequence))

@app.route('/seqResults/<string:sequence>')
def seqResults(sequence):
	return render_template('seqResults.html',sequence=sequence)

@app.route('/coordSearch', methods=['GET','POST'])
def coordSearch():
	coord_form = CoordinateForm(request.form)
	if request.method == 'POST' and coord_form.validate():
		coordinate = coord_form.coordinate.data
		return redirect(url_for('coordResults', coordinate=coordinate))

@app.route('/coordResults/<string:coordinate>')
def coordResults(coordinate):
	return render_template('coordResults.html',coordinate=coordinate)

@app.route('/nameSearch', methods=['GET','POST'])
def nameSearch():
	name_form = GeneNameForm(request.form)
	if request.method == 'POST' and name_form.validate():
		name = name_form.name.data
		return redirect(url_for('nameResults', name=name))

@app.route('/nameResults/<string:name>')
def nameResults(name):
	return render_template('nameResults.html',name=name)

if __name__ == '__main__':
    app.secret_key = 'super_secret_key' # gives us access to the session to flash messages
    app.debug = True
    app.run(host='0.0.0.0', port=5000)