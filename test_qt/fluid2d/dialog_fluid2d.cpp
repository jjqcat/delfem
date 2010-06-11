#include "dialog_fluid2d.h"
#include "ui_dialog_fluid2d.h"

Dialog_Fluid2D::Dialog_Fluid2D(const Fem::Field::CFieldWorld& world_input,
                               Fem::Eqn::CEqnSystem_Fluid2D& fluid_input,
                               QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog_Fluid2D),
    world(world_input),
    fluid(fluid_input)
{
    ui->setupUi(this);

    const unsigned int id_field_velo = fluid.GetIdField_Velo();
    const Fem::Field::CField& field = world.GetField(id_field_velo);
    const std::vector<unsigned int>& aIdEA = field.GetAry_IdElemAry();
    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
    {
        unsigned int id_ea = aIdEA[iiea];
        char str_id[16];
        sprintf(str_id,"%d",id_ea);
        ui->comboBox_IdEA->addItem(str_id);
    }

   connect(ui->comboBox_IdEA,SIGNAL(currentIndexChanged(int)),
           this,SLOT(comboBox_currentIndexChenged(int)));
   connect(ui->doubleSpinBox_myu,SIGNAL(valueChanged(double)),
           this,SLOT(matPorpChanged()));
   connect(ui->doubleSpinBox_rho,SIGNAL(valueChanged(double)),
           this,SLOT(matPorpChanged()));

   this->comboBox_currentIndexChenged(0);
}

Dialog_Fluid2D::~Dialog_Fluid2D()
{
    delete ui;
}

void Dialog_Fluid2D::changeEvent(QEvent *e)
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

void Dialog_Fluid2D::comboBox_currentIndexChenged(int)
{

    unsigned int id_ea0;
    {
        unsigned int ind = ui->comboBox_IdEA->currentIndex();
        QString itemData = ui->comboBox_IdEA->itemText(ind);
        id_ea0 = itemData.toInt();
    }
    const Fem::Eqn::CEqn_Fluid2D& eqn = fluid.GetEqnation(id_ea0);
    ui->doubleSpinBox_myu->setValue(eqn.GetMyu());
    ui->doubleSpinBox_rho->setValue(eqn.GetRho());
    std::cout << "Myu" << " " << eqn.GetMyu() << std::endl;
}

void Dialog_Fluid2D::matPorpChanged()
{
    unsigned int id_ea0;
    {
        unsigned int ind = ui->comboBox_IdEA->currentIndex();
        QString itemData = ui->comboBox_IdEA->itemText(ind);
        id_ea0 = itemData.toInt();
    }
    Fem::Eqn::CEqn_Fluid2D eqn = fluid.GetEqnation(id_ea0);
    eqn.SetMyu(ui->doubleSpinBox_myu->value());
    eqn.SetRho(ui->doubleSpinBox_rho->value());
    fluid.SetEquation(eqn);
}
