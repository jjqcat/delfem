#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "glwidget.h"

#include "dialog_solid2d.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->pushButton_prob,SIGNAL(clicked()),
            this,SLOT(setNewProblem()));
    connect(ui->pushButton_matprop,SIGNAL(clicked()),
            this,SLOT(showDialogMatProp()));

    glWidget = new GLWidget_Solid2d;
    ui->scrollArea->setWidget(glWidget);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::changeEvent(QEvent *e)
{
    QMainWindow::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

void MainWindow::setNewProblem(){
    std::cout << "setNewProblem()" << std::endl;
    glWidget->SetNewProblem();
}

void MainWindow::showDialogMatProp(){
    std::cout << "show Dialog Materila Property()" << std::endl;
    Dialog_Solid2D* dialog = new Dialog_Solid2D(glWidget->world,glWidget->solid,this);
    ((QDialog*)dialog)->setModal(true);
    ((QDialog*)dialog)->exec();
    delete dialog;
}
