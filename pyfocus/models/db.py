# am i optimizing prematurely? do we need this heavy machinary?


from sqlalchemy import Column, Integer, String, Float, ForeignKey, create_engine
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm import relationship, sessionmaker

Base = declarative_base()


def load_db(path):
    # create engine, and ensure that tables exist
    engine = create_engine("sqlite:///{}".format(path))
    Base.metadata.create_all(engine)

    # create a session for all db operations
    session = sessionmaker(bind=engine)

    return session


def build_model(gene_info, snp_info, db_ref_panel, weights, ses, attrs, method):

    mol_feature = MolecularFeature(
        ens_gene_id=gene_info["gene_id"],
        ens_tx_id=gene_info["tx_id"],
        name=gene_info["name"],
        chrom=gene_info["chrom"],
        tx_start=gene_info["tx_start"],
        tx_stop=gene_info["tx_stop"]
    )

    model = Model(
        name="foo",
        inference=method,
        ref_panel=db_ref_panel,
        mol_feature=mol_feature
    )

    # might be slow...
    model.weights = [
        Weight(
            rsid=s_info.ID,
            chrom=s_info.CHR,
            pos=s_info.BP,
            effect_allele=s_info.A1,
            alt_allele=s_info.A0,
            effect=w,
            std_error=se
        ) for w, se, s_info in zip(weights, ses, snp_info.iterrows())
    ]

    model.attrs = [
        ModelAttribute(
            name=name,
            value=value
        ) for name, value in attrs.iteritems()
    ]

    return model


class FocusMixin(object):

    id = Column(Integer, primary_key=True, autoincrement=True)

    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()


class RefPanel(Base, FocusMixin):

    # Main attributes
    name = Column(String(128), nullable=False)
    tissue = Column(String(128), nullable=False)
    assay = Column(String(128))

    # Link to predictive models
    models = relationship("Model")


class Model(Base, FocusMixin):

    # Main attribute
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


class ModelAttribute(Base, FocusMixin):

    # Main
    name = Column(String(128), nullable=False)
    value = Column(Float)

    # Associate with model
    model_id = Column(Integer, ForeignKey("model.id"))
    model = relationship("Model", back_populates="model")


class MolecularFeature(Base, FocusMixin):

    ens_gene_id = Column(String(64), nullable=False)
    ens_tx_id = Column(String(64))

    name = Column(String(64))
    type = Column(String(64))

    chrom = Column(String(10), nullable=False)
    tx_start = Column(Integer, nullable=False)
    tx_stop = Column(Integer, nullable=False)

    models = relationship("Model", back_populates="molecularfeature")


class Weight(Base, FocusMixin):

    rsid = Column(String(128), unique=True, nullable=False)
    chrom = Column(String(2), nullable=False)
    pos = Column(Integer, nullable=False)
    effect_allele = Column(String(32), nullable=False)
    alt_allele = Column(String(32), nullable=False)

    effect = Column(Float, nullable=False)
    std_error = Column(Float)

    model_id = Column(Integer, ForeignKey("model.id"))
    model = relationship("Model", back_populates="model")
