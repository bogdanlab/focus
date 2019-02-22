# am i optimizing prematurely? do we need this heavy machinary?


from sqlalchemy import Column, Integer, String, Float, ForeignKey, create_engine
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm import relationship, sessionmaker

Base = declarative_base()


session = None


def set_session(ssn):
    global session
    session = ssn
    return


def get_session():
    global session
    return session


def load_db(path):
    # create engine, and ensure that tables exist
    engine = create_engine("sqlite:///{}".format(path))
    Base.metadata.create_all(engine)

    # create a session for all db operations
    factory = sessionmaker(bind=engine)

    set_session(factory())
    ssn = get_session()

    return ssn


def build_model(gene_info, snp_info, db_ref_panel, weights, ses, attrs, method):

    mol_feature = MolecularFeature(
        ens_gene_id=gene_info["geneid"],
        ens_tx_id=gene_info["txid"],
        mol_name=gene_info["name"],
        type=gene_info["type"],
        chrom=gene_info["chrom"],
        tx_start=gene_info["txstart"],
        tx_stop=gene_info["txstop"]
    )

    model = Model(
        inference=method,
        ref_panel=db_ref_panel,
        mol_feature=mol_feature
    )

    if ses is None:
        ses = [None] * len(weights)

    # might be slow...
    model.weights = [
        Weight(
            snp=snp_info.iloc[idx].snp,
            chrom=snp_info.iloc[idx].chrom,
            pos=int(snp_info.iloc[idx].pos),
            effect_allele=snp_info.iloc[idx].a1,
            alt_allele=snp_info.iloc[idx].a0,
            effect=w,
            std_error=ses[idx]
        ) for idx, w in enumerate(weights)
    ]

    model.attrs = [
        ModelAttribute(
            attr_name=name,
            value=value
        ) for name, value in attrs.items()
    ]

    return model


class FocusMixin(object):

    id = Column(Integer, primary_key=True, autoincrement=True)

    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()


class RefPanel(Base, FocusMixin):

    # Main attributes
    ref_name = Column(String(128), nullable=False)
    tissue = Column(String(128), nullable=False)
    assay = Column(String(128))

    # Link to predictive models
    models = relationship("Model", back_populates="ref_panel")


class Model(Base, FocusMixin):

    # Main attribute
    inference = Column(String(128), nullable=False)

    # Link back to RefPanel
    ref_id = Column(Integer, ForeignKey('refpanel.id'))
    ref_panel = relationship("RefPanel", back_populates="models")

    # Link to weights for this model
    weights = relationship("Weight", back_populates="model")

    # Link to molecular attributes
    mol_id = Column(Integer, ForeignKey('molecularfeature.id'))
    mol_feature = relationship("MolecularFeature", back_populates="models")

    # Link to general attributes
    attrs = relationship("ModelAttribute", back_populates="model")


class ModelAttribute(Base, FocusMixin):

    # Main
    attr_name = Column(String(128), nullable=False)
    value = Column(Float)

    # Associate with model
    model_id = Column(Integer, ForeignKey("model.id"))
    model = relationship("Model", back_populates="attrs")


class MolecularFeature(Base, FocusMixin):

    ens_gene_id = Column(String(64), nullable=False)
    ens_tx_id = Column(String(64))

    mol_name = Column(String(64))
    type = Column(String(64))

    chrom = Column(String(10), nullable=False)
    tx_start = Column(Integer)
    tx_stop = Column(Integer)

    models = relationship("Model", back_populates="mol_feature")


class Weight(Base, FocusMixin):

    snp = Column(String(128), nullable=False)
    chrom = Column(String(2), nullable=False)
    pos = Column(Integer, nullable=False)
    effect_allele = Column(String(32), nullable=False)
    alt_allele = Column(String(32), nullable=False)

    effect = Column(Float, nullable=False)
    std_error = Column(Float)

    model_id = Column(Integer, ForeignKey("model.id"))
    model = relationship("Model", back_populates="weights")
