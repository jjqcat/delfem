#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class GLWidget_Solid2d;
namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

protected:
    void changeEvent(QEvent *e);
public slots:
    void setNewProblem();
    void showDialogMatProp();

private:
    Ui::MainWindow *ui;
    GLWidget_Solid2d* glWidget;
};

#endif // MAINWINDOW_H
