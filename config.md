### .vimrc
```
set runtimepath+=~/.vim_runtime
let g:go_version_warning = 0
source ~/.vim_runtime/vimrcs/basic.vim
source ~/.vim_runtime/vimrcs/filetypes.vim
source ~/.vim_runtime/vimrcs/plugins_config.vim
source ~/.vim_runtime/vimrcs/extended.vim


inoremap jk <ESC>
set t_Co=256   " This is may or may not needed.
set number
set background=dark
colorscheme PaperColor
set pastetoggle=<F3>
try
source ~/.vim_runtime/my_configs.vim
catch
endtry
```
### Course Links:
1. https://www.mia.uni-saarland.de/mvc/teaching-ws19.shtml
2. https://webmail.uni-saarland.de/
