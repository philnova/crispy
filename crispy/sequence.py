"""
Class to manage interactions with the database
"""
from models import Base, Sequence
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

engine = create_engine('sqlite:///guidesequences.db')
Base.metadata.bind = engine

DBSession = sessionmaker(bind=engine)
session = DBSession()

class NewSequence(object):

	def __init__(self):
		self.sequence = None
		self.assembly = None
		self.organism = None
		self.chromosome = None
		self.start = None
		self.stop = None
		self.direction = None
		self.gene = None
		self.score = None
		self.pam = None

	def __str__(self):
		return self.sequence + self.assembly + self.organism + self.chromosome + str(self.start) + str(self.stop) + self.direction

	def load(self):
		session.add(Sequence(organism = self.organism, assembly = self.assembly, sequence = self.sequence, chromosome = self.chromosome, start = self.start, stop = self.stop, direction = self.direction, score = self.score, gene = self.gene, pam=self.pam))
		session.commit()

