# am i optimizing prematurely? do we need this heavy machinary?

from sqlalchemy import Column, Integer, String, Float, ForeignKey, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

Base = declarative_base()


def load_db(path):
    # TODO: check schema here and make sure everything is copacetic before continuing
    engine = create_engine("sqlite:///{}".format(path))
    session = sessionmaker(bind=engine)

    return session


def append_model(model, session):
    return


class RefPanel(Base):
    __tablename__ = "refpanel"

    # Main attribute
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(128), nullable=False)
    tissue = Column(String(128), nullable=False)
    assay = Column(String(128))

    # Link to predictive models
    models = relationship("Model")


class Model(Base):
    __tablename__ = "model"

    # Main attribute
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(128), nullable=False)
    inference = Column(String(128), nullable=False)

    # Link back to RefPanel
    ref_id = Column(Integer, ForeignKey('refpanel.id'))
    ref_panel = relationship("RefPanel", back_populates="model")

    # Link to weights for this model
    weights = relationship("Weight", back_populates="model")

    # Link to molecular attributes
    mol_id = Column(Integer, ForeignKey('molecularfeature.id'))
    mol_feature = relationship("MolecularFeature", back_populates="model")

    # Link to general attributes
    attrs = relationship("ModelAttribute", back_populates="model")


class ModelAttribute(Base):
    __tablename__ = "modelattr"

    # Main
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(128), nullable=False)
    value = Column(Float)

    # Associate with model
    model_id = Column(Integer, ForeignKey("model.id"))
    model = relationship("Model", back_populates="model")


class MolecularFeature(Base):
    __tablename__ = "molecularfeature"

    id = Column(Integer, primary_key=True, autoincrement=True)

    ens_gene_id = Column(String(64), nullable=False)
    ens_tx_id = Column(String(64))

    name = Column(String(64))
    type = Column(String(64))

    chrom = Column(String(10), nullable=False)
    tx_start = Column(Integer, nullable=False)
    tx_stop = Column(Integer, nullable=False)

    models = relationship("Model", back_populates="molecularfeature")


class Weight(Base):
    __tablename__ = "weight"

    id = Column(Integer, primary_key=True, autoincrement=True)

    rsid = Column(String(128), unique=True, nullable=False)
    chrom = Column(String(10), nullable=False)
    pos = Column(Integer, nullable=False)
    effect_allele = Column(String(32), nullable=False)
    alt_allele = Column(String(32), nullable=False)

    effect = Column(Float, nullable=False)
    std_error = Column(Float)

    model_id = Column(Integer, ForeignKey("model.id"))
    model = relationship("Model", back_populates="model")
