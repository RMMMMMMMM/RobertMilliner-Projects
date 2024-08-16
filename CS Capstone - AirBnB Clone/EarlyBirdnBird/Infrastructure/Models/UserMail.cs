using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations.Schema;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Infrastructure.Models
{
    public class UserMail
    {

        [Key]
        public int Id { get; set; }

        [ForeignKey("UserName")]
        [Display(Name = "From")]
        public AppUser? fromUser { get; set; }

        [Required]
        public string? toUser { get; set; }

        [Required]
        public string? subject { get; set; }

        [Required]
        public string? body { get; set; }

        [Required]
        [DataType(DataType.Date)]
        public DateTime? DateSent { get; set; }

        public bool? isRead { get; set; }

        public bool? isDeleted { get; set;}
    }
}
