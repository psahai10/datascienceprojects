team_summary_df$rush_epa_diff_norm <- team_summary_df$rush_epa_diff*100

pmain <- ggplot(team_summary_df,
                        aes(x = rush_epa_diff, y = Total_Score_Diff)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Rushing EPA/Attempt differential",
       y = "Regular season score differential") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 13))

g <- ggplot(team_summary_df, aes(rush_epa_diff))
g + geom_histogram(binwidth = 0.04)

p <- ggplot(team_summary_df, aes(Total_Score_Diff))
p + geom_histogram(binwidth = 15) +   coord_flip()

xbox <- axis_canvas(pmain, axis = "x") +
  geom_histogram(data=team_summary_df, aes(x =rush_epa_diff), binwidth = 0.03, colour='black', size=0.5)

ybox <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
  geom_histogram(data=team_summary_df, aes(x =Total_Score_Diff), binwidth = 18, colour='black', size=0.5) + coord_flip()

p1 <- insert_xaxis_grob(pmain, xbox, grid::unit(0.8, "in"), position = "top")

p2 <- insert_yaxis_grob(p1, ybox, grid::unit(1, "in"), position = "right")
ggdraw(p2)

library(cowplot)

pmain <- ggplot(data = mpg, aes(x = cty, y = hwy, color = factor(cyl))) + 
  geom_point() + 
  xlab("City driving (miles/gallon)") +
  ylab("Highway driving (miles/gallon)") +
  theme_minimal()

xbox <- axis_canvas(pmain, axis = "x", coord_flip = TRUE) + 
  geom_boxplot(data = mpg, aes(y = cty, x = factor(cyl), color = factor(cyl))) + 
  scale_x_discrete() + coord_flip()

ybox <- axis_canvas(pmain, axis = "y") + 
  geom_boxplot(data = mpg, aes(y = hwy, x = factor(cyl), color = factor(cyl))) +
  scale_x_discrete()

p1 <- insert_xaxis_grob(pmain, xbox, grid::unit(0.6, "in"), position = "top")
p2 <- insert_yaxis_grob(p1, ybox, grid::unit(0.6, "in"), position = "right")
ggdraw(p2)
