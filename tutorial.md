####This tutorial is using [linuxmint-17.3-cinnamon-64bit.iso]()

#### Assume you have a new system.

Make sure you have Java dowloaded or updated. The following command will take a few minutes to complete.

``` $ sudo apt-get install openjdk-7-jre ```

Also important is to install the wrapper for supporting Java applications on your web browser, the Iced-Tea Java Plugin.

``` $ sudo apt-get install icedtea-plugin ```

To make sure you have the most updated version of Java, run the following command:

```$ sudo update-alternatives --config java ```

Optional (for development):
This adds the respository for sublime text 3. If you are using a different
Text editor, then you don't need this repo. Likewize, you don't need to 

```$ sudo add-apt-repository ppa:webupd8team/sublime-text-3```

This updates and upgrades your current system.

```$ sudo apt-get update && sudo apt-get upgrade -y```

Adds the git repository to your files. This is version controll for our project.

```$ sudo apt-get install git -y```

Download the master branch from our prject. Currently it is located at the URI mentioned below.

```$ git clone -b mod_wsgi https://github.com/superphy/semantic.git ~/superphy```

#### Optional:
If you want to use sublime-text as your word processor, now is the time to Download it.

```$ sudo apt-get install sublime-text-installer -y```

- Install the other sudo packages we need to run SuperPhy on the machine

####note: sort these!

```$ sudo apt-get install apache2 curl git libapache2-mod-wsgi libyajl2 MUMmer muscle python-dev python-virtualenv wget xvfb -y```

- To enable mod_wsgi, run the following command:

```$ sudo a2enmod wsgi``` 

- Configure and Enable your virtual host: (deploying to apache)

- copy 000-default.conf to /etc/apache2/sites-available

```$ sudo cp ~/superphy/development_virtualhost.conf /etc/apache2/sites-available/000-default.conf```

```$ sudo service apache2 reload```

***
####Assume you have set up the system, and have just downloading the new

$ cd ~/superphy

$ ln -s ~/superphy ~/Desktop/superphy

start and download the blazegraph java file.

```$ bash superphy/database/scripts/start.sh```

------What's happening here?--------

$ apache='/var/www/SuperPhy'
$ project=$(pwd)/var_www_SuperPhy
$ sudo rm -f $apache
$ sudo ln -s $project $apache

virtualenv --no-site-packages venv

$ source venv/bin/activate
$ pip install --upgrade pip && pip install -r venv/requirements.txt
$ nodeenv -p --prebuilt --requirements=venv/npm-requirements.txt


#### How to deactivate

##### This may become deprecated soon. These are bower js components for mithril.
if ! find var_www_SuperPhy/SuperPhy/static/js/bower_components/mithril-components | read v; then
    git clone --depth=1 https://github.com/eddyystop/mithril-components.git var_www_SuperPhy/SuperPhy/static/js/bower_components/mithril-components;
fi

- This is a very large download. If you aren't going to be uploading data, don't bother downloading it.
```
$ mkdir blast &> /dev/null 
$ cd blast
if ! find ncbi*/ | read v; then
    wget -rc -nd -A "*x64-linux.tar.gz" "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
    tar zxvpf *x64-linux.tar.gz
fi && cd ..
```
#### How to Recompile

bash var_www_SuperPhy/SuperPhy/static/compile.sh
sudo /etc/init.d/apache2 reload
