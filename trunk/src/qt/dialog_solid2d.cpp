#include "dialog_solid2d.h"
#include "ui_dialog_solid2d.h"

Dialog_Solid2D::Dialog_Solid2D(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog_Solid2D)

{
    ui->setupUi(this);
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
