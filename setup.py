import os


if not os.path.exists("scripts.zip"):
    os.system('wget https://github.com/tachengtatangi/unfazed-longreads/raw/main/scripts.zip')
if not os.path.exists("scripts.zip"):
    print("You did not download the scripts zip file!!!\n"
          "please execute this command\n\n"
          "wget https://github.com/tachengtatangi/unfazed-longreads/raw/main/scripts.zip\n\n"
          "to ensure that the file is downloaded")

os.system('unzip scripts.zip')
os.system('rm scripts.zip')
os.system('find ~ -type d -name unfazed | xargs -I {} cp informative_site_finder.py __init__.py __main__.py '
          'read_collector.py site_searcher.py snv_phaser.py sv_phaser.py unfazed.py {}')
os.system('rm informative_site_finder.py __init__.py __main__.py '
          'read_collector.py site_searcher.py snv_phaser.py sv_phaser.py unfazed.py')
print("--------------------------------------\n"
      "Congratulations!\nYou have configured successfully!\n"
      "--------------------------------------")
