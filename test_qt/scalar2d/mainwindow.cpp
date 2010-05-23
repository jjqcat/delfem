#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "glwidget.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->pushButton,SIGNAL(clicked()),
            this,SLOT(setNewProblem()));

    glWidget = new GLWidget;
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
