-- links-to-html.lua
-- https://stackoverflow.com/a/49396058
function Link(el)
  el.target = string.gsub(el.target, "%.md", ".html")
  return el
end

-- -- https://stackoverflow.com/a/48172069
-- -- fix-links-single-file.lua
-- function Link (link)
--   link.target = link.target:gsub('.+%.md%#(.+)', '#%1')
--   return link
-- end
    