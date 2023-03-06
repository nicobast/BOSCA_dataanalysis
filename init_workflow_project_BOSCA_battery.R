require(workflowr)
### see for further instruction: https://jdblischak.github.io/workflowr/articles/wflow-01-getting-started.html

# apply git information for version control
wflow_git_config(user.name = "nicobast", user.email = "nico.bast@kgu.de")

#detect system
ifelse(Sys.info()['sysname']=='Linux',
       home_path<-'~',
       home_path<-'C:/Users/Nico')

wflow_start(paste0(home_path,"/PowerFolders/project_BOSCA_battery/"),existing=T)

#build --> docs ar ebuild and can be viewed on local machine
wflow_build()
wflow_view()


#publish --> can be viewed by others
wflow_publish(c("analysis/index.Rmd", "analysis/about.Rmd", "analysis/license.Rmd",
                "analysis/preprocessing_PDmemory.Rmd","analysis/analysis_PDmemory.Rmd"),
              "Publish the initial files for myproject")

#deploy to github
wflow_use_github("nicobast")
wflow_git_push()

