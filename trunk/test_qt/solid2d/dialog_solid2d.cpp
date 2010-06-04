#include "dialog_solid2d.h"
#include "ui_dialog_solid2d.h"

Dialog_Solid2D::Dialog_Solid2D(const Fem::Field::CFieldWorld& world_input,
                               Fem::Eqn::CEqnSystem_Solid2D& solid_input,
                               QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog_Solid2D),
    world(world_input),
    solid(solid_input)
{
    ui->setupUi(this);

    const unsigned int id_field_disp = solid.GetIdField_Disp();
    const Fem::Field::CField& disp = world.GetField(id_field_disp);
    const std::vector<unsigned int>& aIdEA = disp.GetAry_IdElemAry();
    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
    {
        unsigned int id_ea = aIdEA[iiea];
        char str_id[16];
        sprintf(str_id,"%d",id_ea);
        ui->comboBox_IdEA->addItem(str_id);
    }
    const unsigned int id_ea0 = aIdEA[0];
    const Fem::Eqn::CEqn_Solid2D& eqn = solid.GetEqnation(id_ea0);
    double young, poisson;
    eqn.GetYoungPoisson(young,poisson);
    ui->doubleSpinBox_Young->setValue(young);
    ui->doubleSpinBox_Poisson->setValue(poisson);

   connect(ui->comboBox_IdEA,SIGNAL(currentIndexChanged(int)),
           this,SLOT(comboBox_currentIndexChenged(int)));
   connect(ui->doubleSpinBox_Young,SIGNAL(valueChanged(double)),
           this,SLOT(matPorpChanged()));
}

Dialog_Solid2D::~Dialog_Solid2D()
{
    delete ui;
}

void Dialog_Solid2D::changeEvent(QEvent *e)
{
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}


void Dialog_Solid2D::comboBox_currentIndexChenged(int)
{

    unsigned int id_ea0;
    {
        unsigned int ind = ui->comboBox_IdEA->currentIndex();
        QString itemData = ui->comboBox_IdEA->itemText(ind);
        id_ea0 = itemData.toInt();
    }
    const Fem::Eqn::CEqn_Solid2D& eqn = solid.GetEqnation(id_ea0);
    double young, poisson;
    eqn.GetYoungPoisson(young,poisson);
    ui->doubleSpinBox_Young->setValue(young);
    ui->doubleSpinBox_Poisson->setValue(poisson);
}

void Dialog_Solid2D::matPorpChanged()
{
    unsigned int id_ea0;
    {
        unsigned int ind = ui->comboBox_IdEA->currentIndex();
        QString itemData = ui->comboBox_IdEA->itemText(ind);
        id_ea0 = itemData.toInt();
    }
    Fem::Eqn::CEqn_Solid2D eqn = solid.GetEqnation(id_ea0);
    double young = ui->doubleSpinBox_Young->value();
    double poisson = ui->doubleSpinBox_Poisson->value();
    eqn.SetYoungPoisson(young,poisson,true);
    solid.SetEquation(eqn);
}
