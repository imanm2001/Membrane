#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <QListWidgetItem>
#include <q3d/triangle.h>
#include <q3d/edge.h>
#include <QShowEvent>
#include <QtCore>
#include "Phys/membrane.h"
#include "utils/mthread.h"
QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
public slots:
    void listChanged();
    void pushButtClicked(bool);
    void saveObjButtClicked(bool);
    void captureButtClicked(bool);
    void swapPushButtClicked(bool);
    void randomPushButtClicked(bool);
    void focusCamera(bool);
protected:
    void showEvent(QShowEvent *event);
private:
    bool _init;
    bool _running;
    QMutex  _mutex;
    Ui::MainWindow *ui;
    int _tri1Ind;
    bool updateMembrane();
    Q_INVOKABLE bool updatePlan();

    MThread *_mt,*_mt2;
};
#endif // MAINWINDOW_H
