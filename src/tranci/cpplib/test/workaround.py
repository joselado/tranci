# define a builder object that works as the gtk builder, # but that reads
# from file, it emulates the gtk builder from gtk3

ifile = "input.py" # file from where it reads

class Builder:
  def get_object(self,name):
    class Object:
      """ Fake function"""
      def get_text(self):
        """ Fake function """
        ls = open(ifile,"r").readlines()
        for l in ls: # loop over lines
          if len(l.split())>0:
            if l.split()[0][0]!="#": # if it is not a comment
              if l.split("=")[0].split()[0]==name:
                 out = l.split("=")[1].split()[0] # return number as string
                 if '"' in out: out = out.split('"')[1] 
#                 print out
                 return out # return number as string
#        return 0
      def get_active_text(self):
        return self.get_text()
      def show(self): return None
    return Object()  # return fake object
  def add_from_file(self,dummy):  return None
  def connect_signals(self,dummy): return None


class Gtk:
  def main_quit(self): return None
  def main(self): return None




builder = Builder() # dummy object
gtk = Gtk() # dummy object

