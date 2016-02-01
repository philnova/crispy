## BOILERPLATE: NECESSARY FOR SQLALCHEMY

import sys
from sqlalchemy import Column, ForeignKey, Integer, String, Date, Enum, Float, Table
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine

Base = declarative_base() # object lets SQLAlchemy know that our classes correspond to tables in our database

## DEFINE CLASSES

class Sequence(Base):

	# define tablename
	__tablename__ = 'sequence'

	# define mapping
	id = Column(Integer, primary_key = True)
	organism = Column(String, nullable=False)
	assembly = Column(String, nullable=False)
	sequence = Column(String, nullable = False)
	chromosome = Column(Integer, nullable = False)
	start = Column(Integer, nullable = False)
	stop = Column(Integer, nullable = False)
	direction = Column(Enum('forward','reverse'), nullable=False)
	score = Column(Float)
	

## BOILERPLATE: NECESSARY FOR SQLALCHEMY
#### keep at end of file ####
engine = create_engine('sqlite:///guidesequences.db')

Base.metadata.create_all(engine)