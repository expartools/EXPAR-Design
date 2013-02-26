# -*- coding: utf-8 -*-

"""A simple backend for a TODO app, using Elixir"""

import os
from elixir import *
import time

dbdir=os.path.join(os.path.expanduser("~"),".pyqtodo")
dbfile=os.path.join(dbdir,str(int(time.time()))+"tasks.sqlite")

# It's good policy to have your app use a hidden folder in
# the user's home to store its files. That way, you can
# always find them, and the user knows where everything is.

class Task(Entity):
    """
    A task for your TODO list.
    """

    # By inheriting Entity, we are using Elixir to make this
    # class persistent, Task objects can easily be stored in
    # our database, and you can search for them, change them,
    # delete them, etc.

    using_options(tablename='tasks', allowcoloverride='True')
    # This specifies the table name we will use in the database,
    # I think it's nicer than the automatic names Elixir uses.

    finger_id = Field(Integer, primary_key=True)
    finger_type = Field(Unicode,required=True)
    start = Field(Integer,required=True)
    sequence = Field(Unicode,required=True)
    length = Field(Integer,required=True)
    include = Field(Unicode,default=None,required=False)
    exclude = Field(Unicode,default=None,required=False)
    done = Field(Boolean,default=False,required=False)

    # A task has the following:
    #
    # * A sequence ("Buy groceries"). Always try to use unicode
    #    in your app. Using anything else is *not worth
    #    the trouble*.
    #
    # * A date for when it's due.
    #
    # * A "Done" field. Is it done?
    #
    # * A list of tags. For example, "Buy groceries" could be
    # tagged "Home" and "Important". It's ManyToMany because
    # a task can have many tags and a tag can have many tasks.

    def __repr__(self):
        return "Task: "+self.sequence

    # It's always nicer if objects know how to turn themselves
    # into strings. That way you can help debug your program
    # just by printing them. Here, our groceries task would
    # print as "Task: Buy groceries".

# Since earlier I mentioned Tags, we need to define them too:


saveData=None

# Using a database involves a few chores. I put them
# in the initDB function. Just remember to call it before
# trying to use Tasks or Tags!

class Tritemp(Entity):
    """
    A task for your TODO list.
    """

    # By inheriting Entity, we are using Elixir to make this
    # class persistent, Task objects can easily be stored in
    # our database, and you can search for them, change them,
    # delete them, etc.

    using_options(tablename='tritemp',allowcoloverride='True')
    # This specifies the table name we will use in the database,
    # I think it's nicer than the automatic names Elixir uses.

    pair_id = Field(Integer, primary_key=True)
    finger_id = Field(Integer, required=True)
    type = Field(Unicode,required=True)
    start = Field(Unicode,required=True)
    trigger = Field(Unicode,required=True)
    trig_gen= Field(Unicode,required=True)
    temp = Field(Unicode,required=True)
    tri_length = Field(Integer,required=True)
    temp_bayes_class = Field(Unicode,required=True)
    temp_pwm_class = Field(Unicode,required=True)
    temp_p90_score = Field(Float,required=True)
    temp_diff_score = Field(Float,required=True)
    tri_temp_tm = Field(Float,required=True)
    temp_tm = Field(Float,required=True)
    bonds = Field(Integer,required=True)
    done = Field(Boolean,default=False,required=False)

    # A task has the following:
    #
    # * A sequence ("Buy groceries"). Always try to use unicode
    #    in your app. Using anything else is *not worth
    #    the trouble*.
    #
    # * A date for when it's due.
    #
    # * A "Done" field. Is it done?
    #
    # * A list of tags. For example, "Buy groceries" could be
    # tagged "Home" and "Important". It's ManyToMany because
    # a task can have many tags and a tag can have many tasks.

    def __repr__(self):
        return "Task: "+self.trigger

    # It's always nicer if objects know how to turn themselves
    # into strings. That way you can help debug your program
    # just by printing them. Here, our groceries task would
    # print as "Task: Buy groceries".




def initDB():
    # Make sure ~/.pyqtodo exists
    if not os.path.isdir(dbdir):
        os.mkdir(dbdir)
    # Set up the Elixir internal thingamajigs
    metadata.bind = "sqlite:///%s"%dbfile
    setup_all()
    create_all()

    # This is so Elixir 0.5.x and 0.6.x work
    # Yes, it's kinda ugly, but needed for Debian
    # and Ubuntu and other distros.

    global saveData
    import elixir
    #tarea1=Task(finger_type=u'HTH', start=223122, length=29, sequence=u'TACGAGCCAATAGGCATCCAT', include='waht', exclude='this', done=False)
    print [(t.finger_id,t.finger_type,t.start, t.length, t.sequence,t.include, t.exclude, t.done) for t in Task.query.filter().all()]
    if elixir.__version__ < "0.6":
        saveData=session.flush
    else:
        saveData=session.commit
    saveData()
    return dbfile

# Usually, I add a main() function to all modules that
# does something useful, perhaps run unit tests. In this
# case, it demonstrates our backend's functionality.
# You can try it by running it like this::
#
#   python todo.py

# No detailed comments in this one: study it yourself, it's not complicated!

def main():

    # Initialize database
    initDB()
    # Create two tags
    # Create a few tags and tag them
    tarea1=Task(finger_type=u'HTH', start=223122, length=29, sequence=u'TACGAGCCAATAGGCATCCAT', include=2, exclude=1, done=False)
    saveData()
    print "Tasks with l:"
    print [(t.finger_id,t.finger_type,t.start, t.length, t.sequence,t.include, t.exclude, t.done) for t in Task.query.filter().all()]

if __name__ == "__main__":
    main()
