#ifndef MACOPENEVENT_H
#define MACOPENEVENT_H

#include <QFileOpenEvent>

// for MacOS files are opened throught events
class FileOpenEvent : public QObject
{
  Q_OBJECT
  public:
    FileOpenEvent(QObject * parent = 0) : QObject(parent) {}
    virtual ~FileOpenEvent() {}
  protected:
    bool eventFilter(QObject *obj, QEvent *event)
    {
      if (event->type() == QEvent::FileOpen)
      {
        emit open(static_cast<QFileOpenEvent*>(event)->file());
        return true;
      } else
      {
        return QObject::eventFilter(obj, event);
      }
    }

  signals:
    void open(const QString& file);
};

#endif
